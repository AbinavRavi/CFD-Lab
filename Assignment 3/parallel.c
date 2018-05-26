
#include "parallel.h"

#include<mpi.h>


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void init_parallel(int iproc,
               int jproc,
               int imax,
               int jmax,
               int *myrank,
               int *il,
               int *ir,
               int *jb,
               int *jt,
               int *rank_l,
               int *rank_r,
               int *rank_b,
               int *rankt,
               int *omg_i,
               int *omg_j,
               int num_proc)
{

  *omg_i = (*myrank % iproc) + 1;
  *omg_j = (*myrank % jproc) + 1;

  *il = (*omg_i-1)*(imax/iproc) + 1;
  *ir = (*myrank!=iproc)?((*omg_i)*(imax/iproc)):imax;

  *jb = (*omg_j-1)*(jmax/jproc) + 1;
  *jt = (*myrank!=jproc)?((*omg_j)*(jmax/jproc)):jmax;

  if(*il == 0)
  {
    *rank_l = MPI_PROC_NULL;
  }
  else
  {
    *rank_l = *myrank - 1;
  }

  if(*ir == imax+1)
  {
    *rank_r = MPI_PROC_NULL;
  }
  else
  {
    *rank_r = *myrank + 1;
  }

}
