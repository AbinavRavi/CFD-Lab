
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
{

  *omg_i = (*myrank % iproc) + 1;
  *omg_j = (*myrank % jproc) + 1;

  *il = (*omg_i-1)*(imax/iproc) + 1;
  *ir = (*myrank!=iproc)?((*omg_i)*(imax/iproc)):imax;

  *jb = (*omg_j-1)*(jmax/jproc) + 1;
  *jt = (*myrank!=jproc)?((*omg_j)*(jmax/jproc)):jmax;

  if(*il == 1)
  {
    *l_rank = MPI_PROC_NULL;
  }
  else
  {
    *l_rank = *myrank - 1;
  }

  if(*ir == imax)
  {
    *r_rank = MPI_PROC_NULL;
  }
  else
  {
    *r_rank = *myrank + 1;
  }

  if(*jb == 1)
  {
    *b_rank = MPI_PROC_NULL;
  }
  else
  {
    *b_rank = *myrank - iproc;
  }

  if(*jt == jmax)
  {
    *t_rank = MPI_PROC_NULL;
  }
  else
  {
    *t_rank = *myrank + iproc;
  }

}

void pressure_comm(double **P,int il,int ir,int jb,int jt, int l_rank,int r_rank,int b_rank, int t_rank,double *bufSend, double *bufRecv, MPI_Status *status, int chunk )
{
  //Send to left&  recieve from right
    if(l_rank != MPI_PROC_NULL)
    {
      for (int j = jb; j< (jt-jb);++j)
      {
        bufSend[j-1] = P[il][j];
      }
    }
    MPI_Sendrecv(bufSend, sizeof(double)*(jt-jb), MPI_DOUBLE,l_rank, chunk,bufRecv, sizeof(double)*(jt-jb),               MPI_DOUBLE,r_rank, chunk,MPI_COMM_WORLD, status);
     // Recieve teh pressure value
    if(r_rank != MPI_PROC_NULL)
    {
      for (int j = jb;j< (jt-jb); ++j)
      {
        P[ir+1][j] = bufRecv[j-1];
      }
    }


  // send to right & recieve from left
    if(r_rank != MPI_PROC_NULL)
    {
      for(int j = jb; j< (jt-jb);++j)
      {
        bufSend[j+1] = P[ir][j+1];
      }
    }

    MPI_Sendrecv(bufSend, sizeof(double)*(jt-jb), MPI_DOUBLE,r_rank, chunk,bufRecv, sizeof(double)*(jt-jb),               MPI_DOUBLE,l_rank, chunk,MPI_COMM_WORLD, status);

   if(l_rank != MPI_PROC_NULL)
   {
     for(int j = jb; j<(jt-jb); ++j)
     {
       P[il][j+1] = bufRecv[j+1];
     }
   }
  //send to top recieve from bottom
    if(t_rank != MPI_PROC_NULL)
    {
      for(int i = il; i<(ir-il) ; ++i )
      {
        bufSend[i-1] = P[i][jt+1];
      }
    }
    MPI_Sendrecv(bufSend, sizeof(double)*(jt-jb), MPI_DOUBLE,t_rank, chunk,bufRecv, sizeof(double)*(jt-jb),               MPI_DOUBLE,b_rank, chunk,MPI_COMM_WORLD, status);
    if(b_rank != MPI_PROC_NULL)
    {
      for(int i = ir; i<(ir-il); ++i)
      {
        P[i][jb] = bufRecv[i-1];
      }
    }


  //send to bottom recieve from top
    if(b_rank != MPI_PROC_NULL)
               int *t_rank,
               int *omg_i,
      {
        bufSend[i-1] = P[i][jb];
      }
    }
    MPI_Sendrecv(bufSend, sizeof(double)*(jt-jb), MPI_DOUBLE,b_rank, chunk,bufRecv, sizeof(double)*(jt-jb),               MPI_DOUBLE,t_rank, chunk,MPI_COMM_WORLD, status);
    if(t_rank != MPI_PROC_NULL)
    {
      for(int i = il; i<(ir-il); ++i)
      {
        P[i][jt+1] = bufRecv[i-1];
      }
    }

}

void uv_comm(double **U,
            double **V,
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
            MPI_Status *Status,
            int chunk)

{
  /* Velocity Component U */
  // send from left and recieve from right
  if (l_rank != MPI_PROC_NULL){
        for (int j = jb+1; j <= (jt-jb); ++j)
        {
            bufSend[j-2] = U[il][j];
        }
    }

    MPI_Sendrecv(bufSend, sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, bufRecv, sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);

    if (r_rank != MPI_PROC_NULL){
        for (int j = jb+1; j <=(jt-jb); j++){
             U[(ir-il)+1][j]=bufRecv[j-2];
        }
    }
  // send from right and recieve from left
  if (r_rank != MPI_PROC_NULL){
        for (int j = jb+1; j <= (jt-jb); j++)
        {
            bufSend[j-2] = U[(ir-il)-1][j];
        }
    }
    MPI_Sendrecv(bufSend, sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, bufRecv, sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
    if (l_rank !=MPI_PROC_NULL){
        for (int j= jb+1;j <= (jt-jb);j++){
             U[il-1][j]=bufRecv[j-2];
        }
    }
  // send from top and recieve from bottom
   if (t_rank != MPI_PROC_NULL){
        for (int i = il-1; i <= (ir-il); i++)
        {
            bufSend[i-1] = U[i][jt-jb];
        }
    }
       MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
    if (b_rank !=MPI_PROC_NULL){
        for (int i = il-1; i <= (ir-il); i++){
             U[i][jb]=bufRecv[i-1];
        }
    }
  // send from bottom and recieve from top
    if (b_rank != MPI_PROC_NULL){
        for (int i = il-1; i <= (ir-il); i++)
        {
            bufSend[i-1] = U[i][2];
        }
    }

        MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);

    if (t_rank !=MPI_PROC_NULL){
        for (int i = il-1; i <= (ir-il);i++){
             U[i][(jt-jb)+1]=bufRecv[i-1];
        }
    }

}
