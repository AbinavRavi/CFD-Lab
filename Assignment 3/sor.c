#include "sor.h"
#include <math.h>

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
  MPI_Status status
 )
{
	int i, j;
	double rloc = 0.0;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
	//int chunk = 0;
	int x = ir - il;
	int y = jt - jb;
  int chunk = 0;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //printf("start sor \n");
	/* SOR iteration */
	for(i = 1; i <= x; i++) {
		for(j = 1; j <= y; j++) {
			P[i][j] = (1.0-omg)*P[i][j]
				  + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}

  /* Set left & right global domain boundaries according to Neumann boundary conditions */
	if(l_rank == MPI_PROC_NULL)  /* Only receive/send data from/to right */
	{
		for(j = 1; j <= y; j++) {
			P[0][j] = P[1][j];
		}
	}
	if(r_rank == MPI_PROC_NULL)  /* Only send/receive data to/from left */
	{
		for(j = 1; j <= y; j++) {
			P[x+1][j] = P[x][j];
		}
	}

	/* Set top & bottom global domain boundaries according to Neumann boundary conditions */
	if(t_rank == MPI_PROC_NULL)  /* Only receive/send data from/to bottom */
	{
		for(i = 1; i <= x; i++) {
			P[i][y+1] = P[i][y];
		}
	}
	if(b_rank == MPI_PROC_NULL)  /* Only send/receive data to/from top */
	{
		for(i = 1; i <= x; i++) {
			P[i][0] = P[i][1];
		}
	}
 MPI_Barrier(MPI_COMM_WORLD);
	/* Communicate between processes regarding pressure boundaries */
  pressure_comm(P, il, ir, jb, jt, l_rank, r_rank, b_rank, t_rank, bufSend, bufRecv, &status, chunk);

  MPI_Barrier(MPI_COMM_WORLD);
	/* Compute the residual */
	for(i = 1; i <= x; i++) {
		for(j = 1; j <= y; j++) {
			rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
				  ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}
  *res = sqrt(rloc/(x*y));
}
