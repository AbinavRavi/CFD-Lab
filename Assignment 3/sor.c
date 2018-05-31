#include "sor.h"
#include "parallel.h"
#include <math.h>
#include <mpi.h>

void sor(
  double omg,
  double dx,
  double dy,
   int imax,
  int jmax,
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
  int t_rank
 )
{
	int i, j;
	double rloc = 0.0;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
	double *bufSend = 0;
	double *bufRecv = 0;
	MPI_Status status;
	//int chunk = 0;
	int x = ir - il + 1;
	int y = jt - jb + 1;

	/* Set left & right global domain boundaries according to Neumann boundary conditions */
	if(l_rank == MPI_PROC_NULL && r_rank != MPI_PROC_NULL)  /* Only receive/send data from/to right */
	{
		for(j = 1; j <= y; j++) {
			P[0][j] = P[1][j];
		}
	}
	else if(l_rank != MPI_PROC_NULL && r_rank == MPI_PROC_NULL)  /* Only send/receive data to/from left */
	{
		for(j = 1; j <= y; j++) {
			P[x+1][j] = P[x][j];
		}
	}
	else if(l_rank == MPI_PROC_NULL && r_rank == MPI_PROC_NULL)  /* No bordering processes */
	{
		for(j = 1; j <= y; j++) {
			P[0][j] = P[1][j];
			P[x+1][j] = P[x][j];
		}
	}

	/* Set top & bottom global domain boundaries according to Neumann boundary conditions */
	if(t_rank == MPI_PROC_NULL && b_rank != MPI_PROC_NULL)  /* Only receive/send data from/to bottom */
	{
		for(i = 1; i <= x; i++) {
			P[i][y+1] = P[i][y];
		}
	}
	else if(t_rank != MPI_PROC_NULL && b_rank == MPI_PROC_NULL)  /* Only send/receive data to/from top */
	{
		for(i = 1; i <= x; i++) {
			P[i][0] = P[i][1];
		}
	}
	else if(t_rank == MPI_PROC_NULL && b_rank == MPI_PROC_NULL)  /* No bordering processes */
	{
		for(i = 1; i <= x; i++) {
			P[i][0] = P[i][1];
			P[i][y+1] = P[i][y];
		}
	}

	/* SOR iteration */
	for(i = 1; i <= x; i++) {
		for(j = 1; j <= y; j++) {
			P[i][j] = (1.0-omg)*P[i][j]
				  + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}
 printf("Debug: init_matrix called inside init_uvp P RSS\n");
	/* Communicate between processes regarding pressure boundaries */
	pressure_comm(P, il, ir, jt, jb, l_rank, r_rank, b_rank, t_rank, bufSend, bufRecv, &status);
	


	/* Compute the residual */
	for(i = 1; i <= x; i++) {
		for(j = 1; j <= y; j++) {
			rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
				  ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}

	/* Sum the squares of all local residuals then square root that sum for global residual */
	MPI_Allreduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	*res = sqrt((*res)/(imax*jmax));
}
