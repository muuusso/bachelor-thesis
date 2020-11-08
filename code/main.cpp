#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>

#define LOG(x) cout << x << endl;
using namespace std;

const double PI = 4*atan(1);

const unsigned int SEED = 2;

// totale passi Monte Carlo
const int Mc = 200e3;
// passi di equilibratura
const int E = 10e3;
// passi effettivi calcolo osservabili
const int m = Mc - E;

// numero siti lattice per lato del cubo
const int L = 6;
const int N = L*L*L;

// implementazione operatore modulo
// diverso da operatore remainder % per numeri negativi
int mod (int i, int n)
{
  return ((i % n) + n) % n;
}


// calcolo H della configurazione s
tuple<double,double> computeH(double s[L][L][L], double smc[L][L][L], double J, double h)
{
  double H=0, Hmc=0;

  // divido J per due perchè ji = ij. Conto due volte interazione ij 
  J = J / 2;

  for (int j=0; j < N; j++)
  {
    // interazione destra e sinistra
    H -= J * cos(s[j/L/L%L][j/L%L][j%L] - s[j/L/L%L][j/L%L][mod(j%L-1, L)]);
    H -= J * cos(s[j/L/L%L][j/L%L][j%L] - s[j/L/L%L][j/L%L][mod(j%L+1, L)]);
    // interazione sopra e sotto
    H -= J * cos(s[j/L/L%L][j/L%L][j%L] - s[j/L/L%L][mod(j/L%L-1, L)][j%L]);
    H -= J * cos(s[j/L/L%L][j/L%L][j%L] - s[j/L/L%L][mod(j/L%L+1, L)][j%L]);
    // interazione dentro e fuori
    H -= J * cos(s[j/L/L%L][j/L%L][j%L] - s[mod(j/L/L%L-1, L)][j/L%L][j%L]);
    H -= J * cos(s[j/L/L%L][j/L%L][j%L] - s[mod(j/L/L%L+1, L)][j/L%L][j%L]);

    Hmc -= J * cos(smc[j/L/L%L][j/L%L][j%L] - smc[j/L/L%L][j/L%L][mod(j%L-1, L)]);
    Hmc -= J * cos(smc[j/L/L%L][j/L%L][j%L] - smc[j/L/L%L][j/L%L][mod(j%L+1, L)]);

    Hmc -= J * cos(smc[j/L/L%L][j/L%L][j%L] - smc[j/L/L%L][mod(j/L%L-1, L)][j%L]);
    Hmc -= J * cos(smc[j/L/L%L][j/L%L][j%L] - smc[j/L/L%L][mod(j/L%L+1, L)][j%L]);

    Hmc -= J * cos(smc[j/L/L%L][j/L%L][j%L] - smc[mod(j/L/L%L-1, L)][j/L%L][j%L]);
    Hmc -= J * cos(smc[j/L/L%L][j/L%L][j%L] - smc[mod(j/L/L%L+1, L)][j/L%L][j%L]);

    // interazione campo magnetico esterno
    H -= h * cos(s[j/L/L%L][j/L%L][j%L]); 
    Hmc -= h * cos(smc[j/L/L%L][j/L%L][j%L]); 
  }

  return {H,Hmc};
}


// tuple <double,double> local_molecular_field (double s[L][L][L], int j)
// {
//   double Hj[2] = {0};
  
//   // destra e sinistra                                              
//   Hj[0] -= cos(s[j/L/L%L][j/L%L][mod(j%L-1, L)]);      
//   Hj[0] -= cos(s[j/L/L%L][j/L%L][mod(j%L+1, L)]);      
//   // sopra e sotto                                                  
//   Hj[0] -= cos(s[j/L/L%L][mod(j/L%L-1, L)][j%L]);      
//   Hj[0] -= cos(s[j/L/L%L][mod(j/L%L+1, L)][j%L]);      
//   // dentro e fuori                                                 
//   Hj[0] -= cos(s[mod(j/L/L%L-1, L)][j/L%L][j%L]);      
//   Hj[0] -= cos(s[mod(j/L/L%L+1, L)][j/L%L][j%L]);

//   // destra e sinistra                                              
//   Hj[1] -= sin(s[j/L/L%L][j/L%L][mod(j%L-1, L)]);      
//   Hj[1] -= sin(s[j/L/L%L][j/L%L][mod(j%L+1, L)]);      
//   // sopra e sotto                                                  
//   Hj[1] -= sin(s[j/L/L%L][mod(j/L%L-1, L)][j%L]);      
//   Hj[1] -= sin(s[j/L/L%L][mod(j/L%L+1, L)][j%L]);      
//   // dentro e fuori                                                 
//   Hj[1] -= sin(s[mod(j/L/L%L-1, L)][j/L%L][j%L]);      
//   Hj[1] -= sin(s[mod(j/L/L%L+1, L)][j/L%L][j%L]);

//   return {Hj[0], Hj[1]};
// }


tuple<double,double,double,double,double,double> Metropolis (double T, double h)
{
  double J = 1;

  // vettore configurazione
  // assegno a ogni spin un angolo
  double s[L][L][L];
  // vettore proposta Monte Carlo 
  double smc[L][L][L];

  // grandezza tipico proposta passo mc
  // ottimizzo in fase equilibratura
  double delta = .1 / 180. * PI; // 3°

  // Hamiltoniana, Hamiltoniana proposta mc
  double H, Hmc;

  // variabili decisione passo mc
  double p, eta;
  // passi accettati
  double fw_steps;

  // osservabili
  double Ui[m], Mi[m], Mx, My;
  double U2i[m], M2i[m];

  double Hj[2], sj[2], Hj2;
  int i, j;

  // init configurazione random
  for (i=0; i<N; i++)
  {
    // numero random tra 0 e 2*pi
    s[i/L/L%L][i/L%L][i%L] = 2*PI * ((double)rand() / RAND_MAX);
  }

  // Monte Carlo loop
  for (i=0; i < Mc; i++)
  {
    for (j=0; j<N; j++) 
    {
      // proposta Monte Carlo
      smc[j/L/L%L][j/L%L][j%L] = s[j/L/L%L][j/L%L][j%L];
      smc[j/L/L%L][j/L%L][j%L] += delta*((double) rand()/RAND_MAX - 0.5);
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
      for (j=0; j < N; j++) 
      {
        s[j/L/L%L][j/L%L][j%L] = smc[j/L/L%L][j/L%L][j%L]; 
      }
      // aggiorno valore energia
      H = Hmc;
    }

    // // over-relaxation method nella regione critica
    // if ((T < 2.3) && (T > 2.1) && (h==0) && (L>=4))
    // {
    //   for (j=0; j < N; j=j+12)
    //   {
    //     tie(Hj[0], Hj[1]) = local_molecular_field(s, j);
    //     Hj2 = Hj[0]*Hj[0] + Hj[1]*Hj[1];

    //     sj[0] = cos(s[j/L/L%L][j/L%L][j%L]);
    //     sj[1] = sin(s[j/L/L%L][j/L%L][j%L]);
        
    //     sj[0] = -sj[0] + 2*sj[0]*Hj[0]/Hj2*Hj[0];
    //     sj[1] = -sj[1] + 2*sj[1]*Hj[1]/Hj2*Hj[1];

    //     s[j/L/L%L][j/L%L][j%L] = atan(sj[1]/sj[0]);
    //   }
    // }

    // fase di equilibratura e correzione delta 
    if (i < E)
    {
      if (fw_steps > 500) {fw_steps = 250;}

      // ogni 1000 passi controlla passi accettati su passi fatti
      if ((i > 0) and (i % 500 == 0))
      {
        // cambio delta se passi accettati non in intervallo [35%,65%]
        if ((fw_steps > 325) or (fw_steps < 175))
        {
          delta = delta*(1 + (fw_steps/500. - 0.5));
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

      Mx = 0; My = 0;
      for (j=0; j < N; j++)
      {
        Mx += cos(s[j/L/L%L][j/L%L][j%L]);
        My += sin(s[j/L/L%L][j/L%L][j%L]);
      }

      M2i[i-E] = pow(Mx,2) + pow(My,2);
      Mi[i-E] = sqrt(M2i[i-E]);
    }
  }

  double U=0, U2=0, M=0, M2=0;
  double deltaU=0, deltaU2=0, deltaM=0, deltaM2=0;

  // lunghezza massima calcolo autocorrelazione
  int l = 2000;

  double U0, U20, M0, M20;
  double corrU=0, corrU2=0, corrM=0, corrM2=0;

  for (i=0; i < m; i++)
  {
    // medio valori osservabili
    U += Ui[i] / m; U2 += U2i[i] / m;
    M += Mi[i] / m; M2 += M2i[i] / m;

    // punto iniziale correlazione
    if (i % l == 0)
    {
      U0 = Ui[i]; U20 = U2i[i];
      M0 = Mi[i]; M20 = M2i[i];
    }

    // calcolo correlazione e medio
    corrU += U0 * Ui[i] / m;
    corrU2 += U20 * U2i[i] / m;
    
    corrM += M0 * Mi[i] / m;
    corrM2 += M20 * M2i[i] / m;
  }

  deltaU = sqrt(fabs(corrU - pow(U,2))) / sqrt(m-1);
  deltaU2 = sqrt(fabs(corrU2 - pow(U2,2))) / sqrt(m-1);

  deltaM = sqrt(fabs(corrM - pow(M,2))) / sqrt(m-1);
  deltaM2 = sqrt(fabs(corrM2 - pow(M2,2))) / sqrt(m-1);

  double Cv, deltaCv;
  double Chi, deltaChi;

  Cv = (U2 - U*U) / (N * T*T);
  deltaCv = sqrt(pow(deltaU2,2) + pow(2*U*deltaU,2)) / (N * T*T);

  M = M / N;
  deltaM = deltaM / N; 

  Chi = (M2 - M*M*N*N) / (N * T);
  deltaChi = sqrt(pow(deltaM2*N*N,2) + pow(2*M*deltaM,2)) / (N * T); 

  return {Cv,deltaCv,M,deltaM,Chi,deltaChi};
}


int main ()
{
  srand(SEED);

  ofstream tesi;
  tesi.open("data/" + to_string(L) + ".csv");
  tesi << "Temp,h,Cv,deltaCv,M,deltaM,Chi,deltaChi" << endl;

  double Cv, deltaCv;
  double M, deltaM;
  double Chi, deltaChi;

  double Tc = 2.22;
  double T, h; int i, j;

  // calcolo osservabili al variare di h, a T=Tc 
  for (i = 0; i < 41; i++)
  {
    // campo magnetico esterno costante (rispetto a J)
    h = 0 + i*0.05;

    tie(Cv,deltaCv,M,deltaM,Chi,deltaChi) = Metropolis(Tc, h);

    tesi << Tc  << "," << h        << ",";
    tesi << Cv  << "," << deltaCv  << ",";
    tesi << M   << "," << deltaM   << ",";
    tesi << Chi << "," << deltaChi << endl;
  }

  // calcolo osservabili per h = 0
  for (i = 0; i < 41; i++)
  {
    // temperature per h = 0 in kB*T (rispetto a J)
    T = 1.5 + i*0.05;

    tie(Cv,deltaCv,M,deltaM,Chi,deltaChi) = Metropolis(T, 0);

    tesi << T      << "," << 0.       << ",";
    tesi << Cv     << "," << deltaCv  << ",";
    tesi << M      << "," << deltaM   << ",";
    tesi << Chi    << "," << deltaChi << endl;
  }
  
  // calcolo osservabili a diversi valori di h
  for (j = 0; j < 2; j++)
  {
    h = 1 + j;

    for (i = 0; i < 11; i++)
    {
      // temperatura in kB*T (rispetto a J) 
      T = 1.5 + 0.3*i;

      tie(Cv,deltaCv,M,deltaM,Chi,deltaChi) = Metropolis(T, h);

      tesi << T   << "," << h        << ",";
      tesi << Cv  << "," << deltaCv  << ",";
      tesi << M   << "," << deltaM   << ",";
      tesi << Chi << "," << deltaChi << endl;
    }
  }
}

