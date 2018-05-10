#include <math.h>
#include "uvp.h"
#include"helper.h"
#include"boundary_val.h"
#include <stdio.h>


void calculate_dt(double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V)
{

 double Umax = fabs(U[0][0]);

   for( int i = 1 ; i <= imax ; i++ )
   {
      for( int j = 1 ; j <= jmax ; j++ )
      {
         if ( fabs(U[i][j]) > Umax )
            Umax = fabs(U[i][j]);
      }
   }

  double Vmax = fabs(V[0][0]);

   for( int i = 1 ; i <= imax ; i++ )
   {
      for( int j = 1 ; j <= jmax ; j++ )
      {
         if ( fabs(V[i][j]) > Vmax )
            Vmax = fabs(V[i][j]);
      }
   }

 double x, y, z;

 x = (Re/2.0)*pow(((1.0/pow(dx,2.0))+(1.0/pow(dy,2.0))),-1.0);

 y = dx/(fabs(Umax));

 z = dy/(fabs(Vmax));

 double minA = fmin(x,y);
 double minB = fmin(minA,z);
 if( (tau > 0) && (tau < 1))
 {
   *dt = tau * minB;
 }
  printf("dt calculated \n");
}


void calculate_fg(double Re,
		 double GX, double GY,
		 double alpha,
		 double dt,
		 double dx, double dy,
		 int imax, int jmax,
		 double** U, double** V,
		 double** F, double** G,int **flag)
{

    for(int i=0; i<imax-1; i++){

        for(int j=0; j<jmax; j++){

	if((flag[i][j]&(1<<0))&flag[i+1][j]) //

        F[i][j]=U[i][j]+dt*(
                //Central difference scheme for second derivatives
                (1/Re)*((U[i-1][j]-2*U[i][j]+U[i+1][j])/pow(dx,2.0)+(U[i][j-1]-2*U[i][j]+U[i][j+1])/pow(dy,2.0))
                //Modified Donor Cell method on convective term (d(u^2)/dx)
                -(1/dx)*0.25*(
                        (pow((U[i+1][j]+U[i][j]),2.0) - pow((U[i-1][j]+U[i][j]),2.0))
                        +alpha*(fabs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                             )
                 //Modified Donor Cell method on convective term (d(uv)/dy)
                -(1/dy)*0.25*(
                        ((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])- (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))
                    +alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))
                             )
                //Gravity component in x-direction
                +GX);
	}
    }

    for(int i=0; i<imax; i++){

        for(int j=0; j<jmax-1; j++){

	if((flag[i][j]&(1<<0))&flag[i][j+1])

        G[i][j]=V[i][j]+dt*(
                //Central difference Scheme for second derivatives
                (1/Re)*((V[i-1][j]-2*V[i][j]+ V[i+1][j])/pow(dx,2.0)+(V[i][j-1]-2*V[i][j]+ V[i][j+1])/pow(dy,2.0))
                //Modified Donor cell method for d(uv)/dx
                -(1/dx)*0.25*(
                        ((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])- (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))
                    +alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))
                            )
                //modified Donor cell method for d(v^2)/dy
                -(1/dy)*0.25*(
                        (pow((V[i][j]+V[i][j+1]),2.0) - pow((V[i][j-1]+V[i][j]),2.0))
                        +alpha*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
                            )
                +GY);
        }

    }

 /*set boundary values in case of no-slip/free-slip*/
for(int i = 0; i<imax; ++i)
{
	for(int j = 0; j<jmax; ++j)
	{
		if ( B_O(flag[i][j]) )  F[i][j] = U[i][j];

		if ( B_W(flag[i][j]) )  F[i-1][j] = U[i-1][j];

		if ( B_N(flag[i][j]) )  G[i][j] = V[i][j];

		if ( B_S(flag[i][j]) )  G[i][j-1] = V[i][j-1];

		if ( B_NO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j] = V[i][j]; }

		if ( B_NW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j] = V[i][j]; }

		if ( B_SO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j-1] = V[i][j-1]; }

		if ( B_SW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j-1] = V[i][j-1]; }
		//Outflow
		//if (flag[i][j]&(1<<3) ) P[i][j] = 0;
		//Inflow
		if (flag[i][j]&(1<<4) ) F[i][j] = U[i][j];
	}
}
    printf("F,G calculated \n");

}



void calculate_uv(double dt,double dx,double dy,int imax, int jmax, double**U, double**V,double**F,double**G,double **P,int **flag)
{
	for (int i = 1; i<= imax-1;i++)
	{
		for (int j=1; j<=jmax;j++)
		{
			if((flag[i][j]&(1<<0))&flag[i+1][j])
			//update the U component of velocity
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j]-P[i][j]);
		}
	}

	for (int i = 0; i< imax;i++)
	{
		for (int j=0; j<jmax-1;j++)
		{
			if((flag[i][j]&(1<<0))&flag[i][j+1])
			//update the V component of velocity
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j+1] -P[i][j]);
		}
	}
    printf("U,V calculated \n");

}


void calculate_rs(double dt,
		  double dx,
		  double dy,
		  int imax,
		  int jmax,
		  double **F,
		  double **G,
		  double **RS,int **flag)
{

	for(int i=1; i<imax; i++)
   	{
        	for(int j=1; j<jmax; j++)
		{
			if(flag[i][j]&(1<<0))
			RS[i][j] =  (1/dt)*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy );
		}
	}
    printf("RS calculated \n");
}

void calculate_temp(double **temp,double Pr, double Re, int imax,int jmax,double dx, double dy,double dt, double alpha,double **U,double **V)
{
  double dut_dx = (1/dx)*(((U[i][j])*(temp[i][j]+temp[i+1][j])/2)-((U[i-1][j])*(temp[i-1][j]+temp[i][j])/2))+(alpha/dx)*((fabs(U[i][j]))*((temp[i][j]-temp[i+1][j])/2) - (fabs(U[i-1][j]))*((temp[i-1][j]-temp[i][j])/2);

  double dvt_dy = (1/dy)*(((V[i][j])*(temp[i][j]+temp[i][j+1])/2)-((V[i][j-1])*(temp[i][j-1]+temp[i][j])/2))+(alpha/dy)*((fabs(V[i][j]))*((temp[i][j]-temp[i][j+1])/2) - (fabs(V[i][j-1]))*((temp[i][j-1]-temp[i][j])/2);

  double dt2_dx2 = (temp[i+1][j] - 2*temp[i][j] + temp[i-1][j])/(dx*dx);

  double dt2_dy2 = (temp[i][j+1] - 2*temp[i][j] +temp[i][j-1])/(dy*dy);

  double Z = ((1/Re*Pr)*(dt2_dx2+dt2_dy2) - dut_dx - dvt_dy);

  for (int i = 0; i< imax;i++){
    for (int j=0;j<jmax;j++){
      temp[i][j] = temp[i][j]+ (dt*Z);
    }
  }

  for(int i = 0; i<imax; ++i)
  {
  	for(int j = 0; j<jmax; ++j)
  	{
  		if ( B_O(flag[i][j]) )  temp[i][j] = temp[i+1][j];

  		if ( B_W(flag[i][j]) )  temp[i][j] = temp[i-1][j];

  		if ( B_N(flag[i][j]) )  temp[i][j] = temp[i][j+1];

  		if ( B_S(flag[i][j]) )  temp[i][j] = temp[i][j-1];

  		if ( B_NO(flag[i][j]) ) temp[i][j] = (temp[i][j+1] + temp[i+1][j])/2;

  		if ( B_NW(flag[i][j]) ) temp[i][j] = (temp[i][j+1] + temp[i-1][j])/2;

  		if ( B_SO(flag[i][j]) ) temp[i][j] = (temp[i][j-1] + temp[i+1][j])/2;

  		if ( B_SW(flag[i][j]) ) temp[i][j] = (temp[i][j-1] + temp[i-1][j])/2;

  		if (flag[i][j]&(1<<3) ) temp[i][j] = temp[imax/2][jmax/2];

  		if (flag[i][j]&(1<<4) ) temp[i][j] = 0;
  	}
  }

}
