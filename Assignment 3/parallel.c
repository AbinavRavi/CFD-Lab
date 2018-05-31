
#include "parallel.h"


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
  *omg_j = ((*myrank+1-*omg_i)/iproc)+1;

  *il = (*omg_i-1)*(imax/iproc) + 1;
  *ir = (*omg_i!=iproc)?((*omg_i)*(imax/iproc)):imax;

  *jb = (*omg_j-1)*(jmax/jproc) + 1;
  *jt = (*omg_j!=jproc)?((*omg_j)*(jmax/jproc)):jmax;

  if(*il == 1)      *l_rank = MPI_PROC_NULL;
  else              *l_rank = *myrank - 1;

  if(*ir == imax)   *r_rank = MPI_PROC_NULL;
  else              *r_rank = *myrank + 1;


  if(*jb == 1)      *b_rank = MPI_PROC_NULL;
  else              *b_rank = *myrank - iproc;

  if(*jt == jmax)   *t_rank = MPI_PROC_NULL;
  else              *t_rank = *myrank + iproc;


printf("Thread_id: %d omg_ij: %d%d \nil: %d, ir: %d, jb: %d, jt: %d \n",*myrank,*omg_i,*omg_j, *il,*ir,*jb,*jt);

printf("l_rank: %d, r_rank: %d, b_rank: %d, t_rank: %d \n \n", *l_rank,*r_rank,*b_rank,*t_rank);
}

void pressure_comm(double **P,int il,int ir,int jb,int jt,
                  int l_rank,int r_rank,int b_rank, int t_rank,
              double *bufSend, double *bufRecv, MPI_Status *status, int chunk )
{
  int myrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //Send to left &  recieve from right
  if (l_rank == MPI_PROC_NULL)
  {
    printf("Debug1: %d \n",myrank);
    MPI_Recv(P[ir-il+1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else if (r_rank == MPI_PROC_NULL)
  {
    printf("Debug2: %d \n",myrank);
    MPI_Send(P[1], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else
  {
    printf("Debug3: %d \n",myrank);
    MPI_Sendrecv(P[1], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk,
               P[ir-il+1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk,
               MPI_COMM_WORLD, status);
    chunk++;
  }
  Programm_Sync("Synch-1 \n");
  // send to right & recieve from left

  if (l_rank == MPI_PROC_NULL)
  {
    printf("Debug4: %d \n",myrank);
    MPI_Send(P[ir-il], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else if (r_rank == MPI_PROC_NULL)
  {
    printf("Debug5: %d \n",myrank);
    MPI_Recv(P[0], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else
  {
    printf("Debug6: %d \n",myrank);
    MPI_Sendrecv(P[ir-il], sizeof(double)*(jt-jb), MPI_DOUBLE,r_rank, chunk,
                  P[0], sizeof(double)*(jt-jb),MPI_DOUBLE,l_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }
  Programm_Sync("Synch-2 \n");
  //send to top recieve from bottom
  for (int i = 1; i <=ir-il; ++i)
  {

    bufSend[i] = P[i][jt-jb];
  }

  if (b_rank == MPI_PROC_NULL)
  {
    printf("Debug7: %d \n",myrank);
    MPI_Send(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else if (t_rank == MPI_PROC_NULL)
  {
    printf("Debug8: %d \n",myrank);
    MPI_Recv(bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else
  {
    printf("Debug9: %d \n",myrank);
    MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk,
                  bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }

  for (int i = 1; i <=ir-il; ++i)
  {
    P[i][0] = bufRecv[i];
  }
  Programm_Sync("Synch-3 \n");
  ///send to bottom recieve from top
  for (int i = 1; i <=ir-il; ++i)
  {
    bufSend[i] = P[i][1];
  }

  if (b_rank == MPI_PROC_NULL)
  {
    printf("Debug10: %d \n",myrank);
    MPI_Recv(bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else if (t_rank == MPI_PROC_NULL)
  {
    printf("Debug11: %d \n",myrank);
    MPI_Send(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else
  {
    printf("Debug12: %d \n",myrank);
    MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk,
                bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }

  for (int i = 1; i <=ir-il; ++i)
  {
    P[i][ir-il+1] = bufRecv[i];
  }
  Programm_Sync("Synch-4 \n");
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
            MPI_Status *status,
            int chunk)

{
 int myrank;
 /* Velocity U */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //Send to left &  recieve from right
  if (l_rank == MPI_PROC_NULL)
  {
    printf("Debug1: %d \n",myrank);
    MPI_Recv(U[ir-il+1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else if (r_rank == MPI_PROC_NULL)
  {
    printf("Debug2: %d \n",myrank);
    MPI_Send(U[2], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else
  {
    printf("Debug3: %d \n",myrank);
    MPI_Sendrecv(U[2], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk,
               U[ir-il+1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk,
               MPI_COMM_WORLD, status);
    chunk++;
  }
  Programm_Sync("Synch-1 \n");
  // send to right & recieve from left

  if (l_rank == MPI_PROC_NULL)
  {
    printf("Debug4: %d \n",myrank);
    MPI_Send(U[ir-il-1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else if (r_rank == MPI_PROC_NULL)
  {
    printf("Debug5: %d \n",myrank);
    MPI_Recv(U[0], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else
  {
    printf("Debug6: %d \n",myrank);
    MPI_Sendrecv(U[ir-il-1], sizeof(double)*(jt-jb), MPI_DOUBLE,r_rank, chunk,
                  U[0], sizeof(double)*(jt-jb),MPI_DOUBLE,l_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }
  Programm_Sync("Synch-2 \n");
  //send to top recieve from bottom
  for (int i = 1; i <=ir-il; ++i)
  {

    bufSend[i] = U[i][jt-jb];
  }

  if (b_rank == MPI_PROC_NULL)
  {
    printf("Debug7: %d \n",myrank);
    MPI_Send(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else if (t_rank == MPI_PROC_NULL)
  {
    printf("Debug8: %d \n",myrank);
    MPI_Recv(bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else
  {
    printf("Debug9: %d \n",myrank);
    MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk,
                  bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }

  for (int i = 1; i <=ir-il; ++i)
  {
    U[i][0] = bufRecv[i];
  }
  Programm_Sync("Synch-3 \n");

  ///send to bottom recieve from top
  for (int i = 1; i <=ir-il; ++i)
  {
    bufSend[i] = U[i][1];
  }

  if (b_rank == MPI_PROC_NULL)
  {
    printf("Debug10: %d \n",myrank);
    MPI_Recv(bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else if (t_rank == MPI_PROC_NULL)
  {
    printf("Debug11: %d \n",myrank);
    MPI_Send(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else
  {
    printf("Debug12: %d \n",myrank);
    MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk,
                bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }

  for (int i = 1; i <=ir-il; ++i)
  {
    U[i][ir-il+1] = bufRecv[i];
  }
  Programm_Sync("Synch-4 \n");

  /* Velocity V*/

  int myrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //Send to left &  recieve from right
  if (l_rank == MPI_PROC_NULL)
  {
    printf("Debug1: %d \n",myrank);
    MPI_Recv(V[ir-il+1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else if (r_rank == MPI_PROC_NULL)
  {
    printf("Debug2: %d \n",myrank);
    MPI_Send(V[1], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else
  {
    printf("Debug3: %d \n",myrank);
    MPI_Sendrecv(V[1], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk,
               V[ir-il+1], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk,
               MPI_COMM_WORLD, status);
    chunk++;
  }
  Programm_Sync("Synch-1 \n");
  // send to right & recieve from left

  if (l_rank == MPI_PROC_NULL)
  {
    printf("Debug4: %d \n",myrank);
    MPI_Send(V[ir-il], sizeof(double)*(jt-jb), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else if (r_rank == MPI_PROC_NULL)
  {
    printf("Debug5: %d \n",myrank);
    MPI_Recv(V[0], sizeof(double)*(jt-jb), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else
  {
    printf("Debug6: %d \n",myrank);
    MPI_Sendrecv(V[ir-il], sizeof(double)*(jt-jb), MPI_DOUBLE,r_rank, chunk,
                  P[0], sizeof(double)*(jt-jb),MPI_DOUBLE,l_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }
  Programm_Sync("Synch-2 \n");
  //send to top recieve from bottom
  for (int i = 1; i <=ir-il; ++i)
  {

    bufSend[i] = V[i][jt-jb-1];
  }

  if (b_rank == MPI_PROC_NULL)
  {
    printf("Debug7: %d \n",myrank);
    MPI_Send(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else if (t_rank == MPI_PROC_NULL)
  {
    printf("Debug8: %d \n",myrank);
    MPI_Recv(bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else
  {
    printf("Debug9: %d \n",myrank);
    MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk,
                  bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }

  for (int i = 1; i <=ir-il; ++i)
  {
    V[i][0] = bufRecv[i];
  }
  Programm_Sync("Synch-3 \n");

  ///send to bottom recieve from top

  for (int i = 1; i <=ir-il; ++i)
  {
    bufSend[i] = V[i][2];
  }

  if (b_rank == MPI_PROC_NULL)
  {
    printf("Debug10: %d \n",myrank);
    MPI_Recv(bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);
    chunk++;
  }
  else if (t_rank == MPI_PROC_NULL)
  {
    printf("Debug11: %d \n",myrank);
    MPI_Send(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD);
    chunk++;
  }
  else
  {
    printf("Debug12: %d \n",myrank);
    MPI_Sendrecv(bufSend, sizeof(double)*(ir-il), MPI_DOUBLE, b_rank, chunk,
                bufRecv, sizeof(double)*(ir-il), MPI_DOUBLE, t_rank, chunk,
                  MPI_COMM_WORLD, status);
    chunk++;
  }

  for (int i = 1; i <=ir-il; ++i)
  {
    V[i][ir-il+1] = bufRecv[i];
  }
  Programm_Sync("Synch-4 \n");

}
