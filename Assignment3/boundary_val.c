#include "mpi.h"
#include"boundary_val.h"

void boundaryvalues(int imax, int jmax,double **U,double **V, int l_rank, int r_rank,int b_rank, int t_rank)
{
 if(b_rank == MPI_PROC_NULL)
{
    for (int i =0; i <=imax+1 ; ++i)
	{
        U[i][0]  = - U[i][1];
    }
    for (int i =0; i <=imax ; ++i)
	{
        V[i][0]  = 0;
    }
}

if(t_rank == MPI_PROC_NULL)
{
    for( int i=0;i<= imax+1;++i)
    {
        U[i][jmax+1]  = 2-U[i][jmax];
    }
    for( int i=0;i<= imax;++i)
    {
        V[i][jmax] = 0;
    }
}

if(l_rank == MPI_PROC_NULL)
{
    for(int j=0;j<=jmax;++j)
    {
      U[0][j] = 0;
    }
    for(int j=0;j<=jmax+1;++j)
    {
      V[0][j] = - V[1][j];
    }

}

if(r_rank == MPI_PROC_NULL)
{
    for(int j=0;j<=jmax;++j)
    {
        U[imax+2][j] = 0;
    }
    for(int j=0;j<=jmax+1;++j)
    {
        V[imax+1][j]  = - V[imax][j];
    }
}

}
