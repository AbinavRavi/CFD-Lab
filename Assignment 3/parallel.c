
#include "parallel.h"
#include <mpi.h>


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


#include<mpi.h>

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
               int *rank_t,
               int *omg_i,
               int *omg_j,
               int num_proc)
{


for(int rank = 0;rank<num_proc;rank++)
{
	omg_i[rank] = (rank % iproc) + 1;
	omg_j[rank] = (rank % jproc) + 1;


	il[rank] = (omg_i[rank]-1)*(imax/iproc) + 1;
	ir[rank] = (rank!=iproc)?(omg_i[rank]*(imax/iproc)):imax;
	
	jb[rank] = (omg_j[rank]-1)*(jmax/jproc) + 1;
	jt[rank] = (rank!=jproc)?(omg_j[rank]*(jmax/jproc)):jmax;


 if((il[rank] = = 0)
  {
    rank_l[rank] = MPI_PROC_NULL;
  }
  else
  {
    rank_l[rank] = (rank)-1;
  }
 
/* right boundary */

if((ir[]) = = imax+1)
{
  *rank_r = MPI_PROC_NULL; 
}
else
{
  *rank_r = (*myrank)+1;
}

}


MPI_Init(NULL, NULL);
MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
MPI_comm_rank(MPI_COMM_WORLD, myrank);













  /*Assigning rank to neighbour for each subdomain rank_l,rank_r,rank_b,rank_t*/
    /*left boundary*/
























}
