void calculate_dt(Re,tau,dt,dx,dy,imax,jmax,U,V)
{
min(double a, double b, double c)
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
 
 else if(i<0)
 {
   dt = 0.05;
 }
}
