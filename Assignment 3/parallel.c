
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
               int *l_rank,
               int *r_rank,
               int *b_rank,
               int *t_rank,
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
    *l_rank = MPI_PROC_NULL;
  }
  else
  {
    *l_rank = *myrank - 1;
  }

  if(*ir == imax+1)
  {
    *r_rank = MPI_PROC_NULL;
  }
  else
  {
    *r_rank = *myrank + 1;
  }

  if(*jb == 0)
  {
    *b_rank = MPI_PROC_NULL;
  }
  else
  {
    *b_rank = *myrank + 1;
  }

  if(*jt == jmax+1)
  {
    *t_rank = MPI_PROC_NULL;
  }
  else
  {
    *t_rank = *myrank + 1;
  }

}

void pressure_comm(double **P,int il,int ir,int jb,int jt, int l_rank,int r_rank,int b_rank, int t_rank,double *bufSend, double *bufRecv, MPI_Status *status, int chunk )
{
  //Send to left&  recieve from right
    int sendtag = 0;
    int recvtag = 0;
    int MPI_Sendrecv(P[il], sizeof(double)*(jt-jb), MPI_DOUBLE,
                l_rank, sendtag,
                P[ir], sizeof(double)*(jt-jb), MPI_DOUBLE,
                r_rank, recvtag,
                MPI_COMM_WORLD, status);
    
    sendtag++;
    recvtag++;

  // send to right & recieve from left
    int MPI_Sendrecv(P[ir], sizeof(double)*(jt-jb), MPI_DOUBLE,
                r_rank, sendtag,
                P[il], sizeof(double)*(jt-jb), MPI_DOUBLE,
                l_rank, recvtag,
                MPI_COMM_WORLD, status);
    
    sendtag++;
    recvtag++;
  //send to top recieve from bottom
    int MPI_Sendrecv(P[jt], sizeof(double)*(jt-jb), MPI_DOUBLE,
                t_rank, sendtag,
                P[jb], sizeof(double)*(jt-jb), MPI_DOUBLE,
                b_rank, recvtag,
                MPI_COMM_WORLD, status);
    
    sendtag++;
    recvtag++;

    int MPI_Sendrecv(P[jb], sizeof(double)*(jt-jb), MPI_DOUBLE,
                b_rank, sendtag,
                P[jt], sizeof(double)*(jt-jb), MPI_DOUBLE,
                t_rank, recvtag,
                MPI_COMM_WORLD, status);
    
    sendtag++;
    recvtag++;
  
}