#include <math.h>
double min(double a, double b, double c)
   { 
    if (a < b)
     {
        if (a < c) 
             return a; 
        else 
             return c;
     }

    if (b < c)
        return b;
    else return c;
   }

void calculate_dt(double Re,
                  double tau,
                  double &dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V)
{

 double Umax = U[0][0];
 
   for( i = 1 ; i < imax ; i++ )
   {
      for( j = 1 ; j < jmax ; j++ )
      {
         if ( U[i][j] > Umax )
            Umax = U[i][j];
      }
   }
   
  double Vmax = V[0][0];
 
   for( i = 1 ; i < imax ; i++ )
   {
      for( j = 1 ; j < jmax ; j++ )
      {
         if ( V[i][j] > Vmax )
            Vmax = V[i][j];
      }
   }
 
 double x, y, z;
 
 x = (Re/2)*(((1/((dx)^2))+(1/((dy)^2)))^(-1));
 y = dx/(fabs(Umax));
 z = dy/(fabs(Vmax));
 
 if( (tau > 0) && (tau < 1))
 {
   dt = tau * min(x,y,z);
 }
 
 
}
