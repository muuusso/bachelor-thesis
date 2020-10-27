/* 
Take angle as input
return abs(sin(angle)) as output
*/

#include<stdio.h>
#include<math.h> /* has  sin(), abs(), and fabs() */

int main(void)
{ 
  for (int i; i < 10; i++)
  {
    printf("Sin(%f) is %f \n", i/10.0, sin(i/10.0));
    printf("Cos(%f) is %f \n", i/10.0, cos(i/10.0));
  }
}
