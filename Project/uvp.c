#include <math.h>
#include "uvp.h"
#include"helper.h"
#include"boundary_val.h"
#include <stdio.h>

/////////////////////////////////////////////////////////////////////
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

   for( int i = 0 ; i < imax ; i++ )
   {
      for( int j = 0 ; j < jmax ; j++ )
      {
         if ( fabs(U[i][j]) > Umax )
            Umax = fabs(U[i][j]);
      }
   }

  double Vmax = fabs(V[0][0]);

   for( int i = 0 ; i < imax ; i++ )
   {
      for( int j = 0 ; j < jmax ; j++ )
      {
         if ( fabs(V[i][j]) > Vmax )
            Vmax = fabs(V[i][j]);
      }
   }

 double x, y, z;

 x = (Re/2.0)*pow(((1.0/pow(dx,2.0))+(1.0/pow(dy,2.0))),-1.0);

 y = dx/(fabs(Umax));

 z = dy/(fabs(Vmax));

  double min = fmin(x,y);
    min = fmin(min,z);

 if( (tau > 0) && (tau < 1))
 {
   *dt = tau * min;
 }
}

/////////////////////////////////////////////////////////////////////////////////
void calculate_fg(double Re,
		 double GX, double GY,
		 double alpha,
		 double dt,
		 double dx, double dy,
		 int imax, int jmax,
		 double** U, double** V,
		 double** F, double** G,int **flag)
{



 /*set boundary values in case of no-slip/free-slip/and inflow*/
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

    for(int i=0; i<imax-1; i++)
    {
        for(int j=0; j<jmax; j++)
		{
		if( ((flag[i][j]&(1<<0))&flag[i+1][j]) || ( (flag[i+1][j] & (1<<3)) && (flag[i][j]&(1<<0))) )
		{
			
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
                             )+GX);

		}
		}
    }


    for(int i=0; i<imax; i++)
	{
        for(int j=0; j<jmax-1; j++)
		{
		if((flag[i][j]&(1<<0))&flag[i][j+1])
		{
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
                            )+GY);

		}
        }

    }

}


////////////////////////////////////////////////////////////////////////////////
void calculate_uv(double dt,double dx,double dy,int imax, int jmax,
		 double**U, double**V,double**F,double**G,double **P,int **flag)
{
	for (int i = 0; i< imax-1;i++)
	{
		for (int j = 0; j<jmax;j++)
		{
			if(((flag[i][j]&(1<<0))&flag[i+1][j]) || ( (flag[i+1][j] & (1<<3)) && (flag[i][j]&(1<<0))))
			//update the U component of velocity
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j]-P[i][j]);
		}
	}

	for (int i = 0; i< imax;i++)
	{
		for (int j = 0; j<jmax-1;j++)
		{
			if((flag[i][j]&(1<<0))&flag[i][j+1])
			//update the V component of velocity
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j+1] -P[i][j]);
		}
	}
    //printf("U,V calculated \n");

}

/////////////////////////////////////////////////////////////////
void calculate_rs(double dt,
		  double dx,
		  double dy,
		  int imax,
		  int jmax,
		  double **F,
		  double **G,
		  double **RS,int **flag)
{

	for(int i=0; i<imax; i++)
   	{
        	for(int j=0; j<jmax; j++)
		{
			if(flag[i][j]&(1<<9))
			RS[i][j] =  (1/dt)*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy );
		}
	}
    //printf("RS calculated \n");
}


void nullify_obstacles(double **U, double **V, double **P, int **flag, int imax, int jmax)
{

	for (int i = 0; i< imax;i++)
	{
	    for (int j=0;j<jmax; j++)
		{
			if(flag[i][j]&((1<<1)|(1<<2)|(1<<11)) )
			{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;
			}

		}
	}

}

void set_gravity(double *gx, double *gy, double t, int prob){

	switch (prob)
	{
	case 1:
		//sudden break;
		if (t<3.0) {
			*gx = 0.75;
		}
		else
		{
			*gx = 0;
		}

		//uphill, downhill

		//speed breaker

		//banking

		//uneven roads

		break;

	case 2:
		*gx = 1;
		break;

	default:
		break;
	}

}


double Force_x(int imax, int jmax, double dely, double **P, int **flag){

	double fx = 0;
	for(int i = 1; i <imax-1; i++)
	{
		for(int j =1; j <jmax-1 ; j++)
		{
			if(flag[i+1][j]&(1<<7))		fx += P[i][j]*dely;
			
			if(flag[i-1][j]&(1<<8))		fx -= P[i][j]*dely;

		}
	}

	return fx;
}


double Force_y(int imax, int jmax, double delx, double **P, int **flag){

	double fy = 0;
	for(int i = 1; i <imax-1; i++)
	{
		for(int j =1; j <jmax-1 ; j++)
		{
			if(flag[i][j+1]&(1<<6))		fy += P[i][j]*delx;
			
			if(flag[i][j-1]&(1<<5))		fy -= P[i][j]*delx;

		}
	}

	return fy;
}

double KE(int imax, int jmax, double **U, double **V, int **flag){

	double ke = 0;
	for(int i = 1; i <imax-1; i++)
	{
		for(int j =1; j <jmax-1 ; j++)
		{
			if(flag[i][j]&(1<<0))		
			
			ke += 0.5*(U[i][j]*U[i][j] + V[i][j]*V[i][j]);

		}
	}

	return ke;
}