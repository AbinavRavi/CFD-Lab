#include <math.h>
#include "uvp.h"
#include"helper.h"
#include "mpi.h"

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

}


void calculate_fg(double Re,
		 double GX, double GY,
		 double alpha,
		 double dt,
		 double dx, double dy,
		 int imax, int jmax,
		 double** U, double** V,
		 double** F, double** G,
     int b_rank, int t_rank, int l_rank, int r_rank
   )
{

    for(int i=1; i<=imax-1; i++){

        for(int j=1; j<=jmax; j++){

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

	for(int i=1; i<=imax; i++){

        for(int j=1; j<=jmax-1; j++){
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

  //Boundary values
  if(b_rank == MPI_PROC_NULL)
 {
      for (int i =1; i <=imax ; ++i)
        G[i][0] = V[i][0];
 }

 if(t_rank == MPI_PROC_NULL)
 {
     for( int i=1;i<= imax;++i)
       G[i][jmax]=V[i][jmax];
 }

 if(l_rank == MPI_PROC_NULL)
 {
     for(int j=1;j<=jmax;++j)
        F[0][j]=U[0][j];
 }

 if(r_rank == MPI_PROC_NULL)
 {
     for(int j=1;j<=jmax;++j)
      F[imax][j]=U[imax][j];
 }

}

void calculate_uv(double dt,double dx,double dy,int imax, int jmax, double**U, double**V,double**F,double**G,double **P)
{
	for (int i = 1; i<= imax-1;i++)
	{
		for (int j=1; j<=jmax;j++)
		{
			//update the U component of velocity
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j]-P[i][j]);
		}
	}

	for (int i = 1; i<= imax;i++)
	{
		for (int j=1; j<=jmax-1;j++)
		{
			//update the V component of velocity
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j+1] -P[i][j]);
		}
	}
}


void calculate_rs(double dt,
		  double dx,
		  double dy,
		  int imax,
		  int jmax,
		  double **F,
		  double **G,
		  double **RS)
{

for(int i=1; i<=imax; i++)
   {

        for(int j=1; j<=jmax; j++)
		{

		RS[i][j] =  (1/dt)*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy );

		}

	}

}
