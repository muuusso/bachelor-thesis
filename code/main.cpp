#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

#define LOG(x) cout << x << endl;
using namespace std;

const double PI = 4*atan(1);

// numero particelle
const int N = 512;
// siti per lato del lattice 3-dimensionale
const int n = cbrt(N);

// totale passi Monte Carlo
const int M = 128e3;
// passi di equilibratura
const int E = 32e3;
// passi effettivi calcolo osservabili
const int m = M - E;


// implementazione operatore modulo
// diverso da operatore remainder % per numeri negativi
int mod (int i, int n)
{
  return ((i % n) + n) % n;
}


// calcolo H della configurazione s
tuple<double,double> computeH(double s[n][n][n], double smc[n][n][n], double J, double h)
{
  double H=0, Hmc=0;

  // divido J per due perchè ji = ij. Conto due volte interazione ij 
  J = J / 2;

  for (int j=0; j < N; j++)
  {
    // interazione destra e sinistra
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][j/n%n][mod(j%n-1, n)]);
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][j/n%n][mod(j%n+1, n)]);
    // interazione sopra e sotto
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][mod(j/n%n-1, n)][j%n]);
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[j/n/n%n][mod(j/n%n+1, n)][j%n]);
    // interazione dentro e fuori
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[mod(j/n/n%n-1, n)][j/n%n][j%n]);
    H -= J * cos(s[j/n/n%n][j/n%n][j%n] - s[mod(j/n/n%n+1, n)][j/n%n][j%n]);

    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][j/n%n][mod(j%n-1, n)]);
    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][j/n%n][mod(j%n+1, n)]);

    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][mod(j/n%n-1, n)][j%n]);
    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[j/n/n%n][mod(j/n%n+1, n)][j%n]);

    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[mod(j/n/n%n-1, n)][j/n%n][j%n]);
    Hmc -= J * cos(smc[j/n/n%n][j/n%n][j%n] - smc[mod(j/n/n%n+1, n)][j/n%n][j%n]);

    // interazione campo magnetico esterno
    H -= h * cos(s[j/n/n%n][j/n%n][j%n]); 
    Hmc -= h * cos(smc[j/n/n%n][j/n%n][j%n]); 
  }

  return {H,Hmc};
}


tuple<double,double,double,double,double,double> Metropolis (double T, double h)
{
  double J = 1;

  // vettore configurazione
  // assegno a ogni spin un angolo
  double s[n][n][n];
  // vettore proposta Monte Carlo 
  double smc[n][n][n];

  // grandezza tipico proposta passo mc
  // ottimizzo in fase equilibratura
  double delta = 3. / 180. * PI; // 3°

  // Hamiltoniana, Hamiltoniana proposta mc
  double H, Hmc;

  // variabili decisione passo mc
  double p, eta;
  // passi accettati
  double fw_steps;

  // osservabili
  double Ui[m], Mi[m][2];
  double U2i[m], M2i[m];

  // init configurazione random
  for (int i=0; i<N; i++)
  {
    // numero random tra 0 e 2*pi
    s[i/n/n%n][i/n%n][i%n] = 2*PI * ((double)rand() / RAND_MAX);
  }

  // Monte Carlo loop
  for (int i=0; i < M; i++)
  {
    for (int j=0; j < N; j++)
    {
      // proposta Monte Carlo
      smc[j/n/n%n][j/n%n][j%n] = s[j/n/n%n][j/n%n][j%n];
      smc[j/n/n%n][j/n%n][j%n] += delta*((double) rand()/RAND_MAX - 0.5);
    }

    tie(H, Hmc) = computeH(s, smc, J, h);

    // estraggo numero random tra 0 e 1
    eta = (double)rand() / RAND_MAX;
    // calcolo "probabilità" transizione
    p = exp(-(Hmc - H) / T); 

    // controllo passo
    if (p > eta)
    {
      // passo accettato, aumento counter
      fw_steps++;
      // aggiorno configurazione
      for (int j=0; j < N; j++) 
      {
        s[j/n/n%n][j/n%n][j%n] = smc[j/n/n%n][j/n%n][j%n]; 
      }
      // aggiorno valore energia
      H = Hmc;
    }

    // fase di equilibratura e correzione delta 
    if (i < E)
    {
      if (fw_steps > 1000) {fw_steps = 500;}

      // ogni 1000 passi controlla passi accettati su passi fatti
      if ((i > 0) and (i % 1000 == 0))
      {
        // cambio delta se passi accettati non in intervallo [35%,65%]
        if ((fw_steps > 650) or (fw_steps < 350))
        {
          delta = delta*(1 + 0.25*(fw_steps/1000. - 0.5));
        }
        // reset counter
        fw_steps = 0;
      }
    }

    // equilibrio
    if (i >= E)
    {
      // calcolo valori osservabili passo corrente (i)
      Ui[i-E] = H; U2i[i-E] = H*H;

      Mi[i-E][0] = 0; Mi[i-E][1] = 0;
      for (int j=0; j < N; j++)
      {
        Mi[i-E][0] += cos(s[j/n/n%n][j/n%n][j%n]);
        Mi[i-E][1] += sin(s[j/n/n%n][j/n%n][j%n]);
      }

      M2i[i-E] = pow(Mi[i-E][0],2) + pow(Mi[i-E][1],2);
    }
  }

  double U=0, U2=0, Mx=0, My=0, M2=0;
  double deltaU=0, deltaU2=0, deltaMx=0, deltaMy=0, deltaM2=0;

  // lunghezza massima calcolo autocorrelazione
  int l = 2000;

  double U0, U20, Mx0, My0, M20;
  double corrU=0, corrU2=0, corrMx=0, corrMy=0, corrM2=0;

  for (int i=0; i < m; i++)
  {
    // medio valori osservabili
    U += Ui[i] / m; U2 += U2i[i] / m;
    Mx += Mi[i][0] / m;
    My += Mi[i][1] / m;
    M2 += M2i[i] / m;

    // punto iniziale correlazione
    if (i % l == 0)
    {
      U0 = Ui[i]; U20 = U2i[i];
      Mx0 = Mi[i][0];
      My0 = Mi[i][1];
      M20 = M2i[i];
    }

    // calcolo correlazione e medio
    corrU += U0 * Ui[i] / m;
    corrU2 += U20 * U2i[i] / m;
    
    corrMx += Mx0 * Mi[i][0] / m;
    corrMy += My0 * Mi[i][1] / m;
    corrM2 += M20 * M2i[i] / m;
  }

  deltaU = sqrt(fabs(corrU - pow(U,2))) / sqrt(m-1);
  deltaU2 = sqrt(fabs(corrU2 - pow(U2,2))) / sqrt(m-1);

  deltaMx = sqrt(fabs(corrMx - pow(Mx,2))) / sqrt(m-1);
  deltaMy = sqrt(fabs(corrMy - pow(My,2))) / sqrt(m-1);
  deltaM2 = sqrt(fabs(corrM2 - pow(M2,2))) / sqrt(m-1);

  double Cv, deltaCv;
  double M, deltaM;
  double Chi, deltaChi;

  Cv = (U2 - U*U) / (N * T*T);
  deltaCv = sqrt(pow(deltaU2,2) + pow(2*U*deltaU,2)) / (N * T*T);

  M = sqrt(Mx*Mx + My*My) / N;
  deltaM = sqrt(Mx*Mx*deltaMx*deltaMx + My*My*deltaMy*deltaMy) / (N*N * M);

  Chi = (M2 - M*M) / (N * T);
  deltaChi = sqrt(pow(deltaM2,2) + pow(2*M*deltaM,2)) / (N*T); 

  return {Cv,deltaCv,M,deltaM,Chi,deltaChi};
}


int main ()
{
  ofstream tesi;
  tesi.open("tesi.csv");
  tesi << "Temp,h,Cv,deltaCv,M,deltaM,Chi,deltaChi" << endl;

  double Cv, deltaCv;
  double M, deltaM;
  double Chi, deltaChi;

  // campo magnetico esterno costante (rispetto a J)
  // range di h / J; J = 1
  double h[3] = {0, 1, 2};

  // temperatura in kB*T (rispetto a J) 
  // range di kB*T/J; J = 1
  double T[40] = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 
                  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3,
                  3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4,
                  4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5};

  for (int j=0; j < 3; j++)
  {
    for (int i=0; i < 40; i++)
    {
      tie(Cv,deltaCv,M,deltaM,Chi,deltaChi) = Metropolis(T[i], h[j]);

      tesi << T[i]   << "," << h[j]     << ",";
      tesi << Cv     << "," << deltaCv  << ",";
      tesi << M      << "," << deltaM   << ";";
      tesi << Chi    << "," << deltaChi << endl;
    }
  }
}
