#include <math.h>
#include "uvp.h"
#include"helper.h"
#include"mpi.h"


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
 double Umax1;
   for( int i = 0 ; i <= imax+2 ; i++ )
   {
      for( int j = 0 ; j <= jmax+1 ; j++ )
      {
         if ( fabs(U[i][j]) > Umax )
            Umax = fabs(U[i][j]);
      }
   }
   MPI_Allreduce(&Umax, &Umax1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  double Vmax = fabs(V[0][0]);
  double Vmax1;
   for( int i = 0 ; i <= imax+1 ; i++ )
   {
      for( int j = 0 ; j <= jmax+2 ; j++ )
      {
         if ( fabs(V[i][j]) > Vmax )
            Vmax = fabs(V[i][j]);
      }
   }
MPI_Allreduce(&Vmax, &Vmax1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

 double x, y, z;

 x = (Re*0.5)*pow( (1.0/(dx*dx)+ 1.0/(dy*dy)),-1.0);

 y = dx/(fabs(Umax1));

 z = dy/(fabs(Vmax1));

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
     int l_rank, int r_rank, int b_rank, int t_rank)
{

for (int i = 0; i <=imax+1; i++) {
        //Boundary conditions for G
        G[i][0] = V[i][0];
        G[i][1] = V[i][1];
	      G[i][jmax+1] = V[i][jmax+1];
    }

    //Boundary conditions for F
    for (int j = 0; j <=jmax+1; j++) {
        F[0][j]=U[0][j];
        F[1][j] = U[1][j];
        F[imax+1][j] = U[imax+1][j];
    }

    for(int i=2; i<=imax; i++){

        for(int j=1; j<=jmax; j++){

        F[i][j]=U[i][j]+dt*(
                //Central difference scheme for second derivatives
                (1/Re)*((U[i-1][j]-2*U[i][j]+U[i+1][j])/(dx*dx)+(U[i][j-1]-2*U[i][j]+U[i][j+1])/(dy*dy))
                //Modified Donor Cell method on convective term (d(u^2)/dx)
                -(1/dx)*0.25*(
                        ( (U[i+1][j]+U[i][j])*(U[i+1][j]+U[i][j]) - (U[i-1][j]+U[i][j])*(U[i-1][j]+U[i][j]) )
                        +alpha*(fabs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                             )
                 //Modified Donor Cell method on convective term (d(uv)/dy)
                -(1/dy)*0.25*(
                        ((V[i-1][j+1]+V[i][j+1])*(U[i][j]+U[i][j+1])- (V[i-1][j]+V[i][j])*(U[i][j-1]+U[i][j]))
                    +alpha*(fabs(V[i-1][j+1]+V[i][j+1])*(U[i][j]-U[i][j+1])-fabs(V[i-1][j]+V[i][j])*(U[i][j-1]-U[i][j]))
                             )
                //Gravity component in x-direction
                +GX);
		}
	}

	for(int i=1; i<=imax; i++){

        for(int j=2; j<=jmax; j++){
        G[i][j]=V[i][j]+dt*(
                //Central difference Scheme for second derivatives
                (1/Re)*((V[i-1][j]-2*V[i][j]+ V[i+1][j])/(dx*dx)+(V[i][j-1]-2*V[i][j]+ V[i][j+1])/(dy*dy))
                //Modified Donor cell method for d(uv)/dx
                -(1/dx)*0.25*(
                        ((U[i+1][j-1]+U[i+1][j])*(V[i][j]+V[i+1][j])- (U[i][j-1]+U[i][j])*(V[i-1][j]+V[i][j]))
                    +alpha*(fabs(U[i+1][j-1]+U[i+1][j])*(V[i][j]-V[i+1][j])-fabs(U[i][j-1]+U[i][j])*(V[i-1][j]-V[i][j]))
                            )
                //modified Donor cell method for d(v^2)/dy
                -(1/dy)*0.25*(
                        ( (V[i][j]+V[i][j+1])*(V[i][j]+V[i][j+1]) - (V[i][j-1]+V[i][j])*(V[i][j-1]+V[i][j]) )
                        +alpha*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
                            )
                +GY);
        }
	}


}



void calculate_uv(double dt,double dx,double dy,int imax, int jmax, double**U, double**V,double**F,double**G,double **P)
{
	for (int i = 1; i<= imax+1;i++)
	{
		for (int j=1; j<=jmax;j++)
		{
			//update the U component of velocity
			U[i][j] = F[i][j] - (dt/dx)*(P[i][j]-P[i-1][j]);
		}
	}

	for (int i = 1; i<= imax;i++)
	{
		for (int j=1; j<=jmax+1;j++)
		{
			//update the V component of velocity
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j] -P[i][j-1]);
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
		RS[i][j] =  (1/dt)*( (F[i+1][j]-F[i][j])/dx + (G[i][j+1]-G[i][j])/dy );
		}

	}

}
