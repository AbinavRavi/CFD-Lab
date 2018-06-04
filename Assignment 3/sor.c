#include "sor.h"
#include <math.h>
#include "parallel.h"

void sor(
  double omg,
  double dx,
  double dy,
  double **P,
  double **RS,
  double *res,
  int il,
  int ir,
  int jb,
  int jt,
  int l_rank,
  int r_rank,
  int b_rank,
  int t_rank,
  double *bufSend,
  double *bufRecv,
  int imax,
  int jmax,
  MPI_Status status
 )
{
	int i, j;
	double rloc = 0.0;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
	int xdim = ir - il+1;
	int ydim = jt - jb+1;
 	int chunk = 0;
  
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	/* Compute the residual */


		/* SOR iteration */
	for(i = 1; i <= xdim; i++) {
		for(j = 1; j <= ydim; j++) {
			P[i][j] = (1.0-omg)*P[i][j]
				  + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}
/* Communicate between processes regarding pressure boundaries */
  pressure_comm(P, il, ir, jb, jt, l_rank, r_rank, b_rank, t_rank, bufSend, bufRecv, &status, chunk);
	
	/* Set left & right global domain boundaries according to Neumann boundary conditions */
	if(l_rank == MPI_PROC_NULL)
	{
		for(j = 1; j <= ydim; j++) {
			P[0][j] = P[1][j];
		}
	}
	if(r_rank == MPI_PROC_NULL)
	{
		for(j = 1; j <= ydim; j++) {
			P[xdim+1][j] = P[xdim][j];
		}
	}

	/* Set top & bottom global domain boundaries according to Neumann boundary conditions */
	if(t_rank == MPI_PROC_NULL)
	{
		for(i = 1; i <= xdim; i++) {
			P[i][ydim+1] = P[i][ydim];
		}
	}
	if(b_rank == MPI_PROC_NULL)
	{
		for(i = 1; i <= xdim; i++) {
			P[i][0] = P[i][1];
		}
	}


	for(i = 1; i <= xdim; i++) {
		for(j = 1; j <= ydim; j++) {
			rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
				  ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}
	//printf("res:%d = %f \n",myrank,rloc);
  MPI_Allreduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *res = sqrt(*res/(imax*jmax));
}
