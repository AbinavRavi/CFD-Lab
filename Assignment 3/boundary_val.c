#include"mpi.h"
#include"boundary_val.h"

void boundaryvalues(int imax, int jmax,double **U,double **V,int b_rank, int t_rank, int l_rank, int r_rank)
{
 if(b_rank == MPI_PROC_NULL)
{
        for (int i =1; i <=imax ; ++i)
	{
        U[i][0]       = - U[i][1];
        V[i][1]       = 0;
    }
}

if(t_rank == MPI_PROC_NULL)
{
    for( int i=1;i<= imax;++i)
    {
        U[i][jmax+1]  = 2-U[i][jmax];
        V[i][jmax]    = 0;
    }
}

if(l_rank == MPI_PROC_NULL)
{
    for(int j=1;j<=jmax;++j)
    {
      U[0][j] = 0;
      V[0][j] = - V[1][j];
    }
}

if(r_rank == MPI_PROC_NULL)
{
    for(int j=1;j<=jmax;++j)
    {
        U[imax][j] = 0;
        V[imax+1][j]  = - V[imax][j];
    }
}

}
