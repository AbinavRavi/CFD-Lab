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


MPI_Barrier(MPI_COMM_WORLD);
printf("Thread_id: %d omg_ij: %d%d \nil: %d, ir: %d, jb: %d, jt: %d \n",*myrank,*omg_i,*omg_j, *il,*ir,*jb,*jt);

printf("l_rank: %d, r_rank: %d, b_rank: %d, t_rank: %d \n \n", *l_rank,*r_rank,*b_rank,*t_rank);
}

void pressure_comm(double **P,int il,int ir,int jb,int jt,
                  int l_rank,int r_rank,int b_rank, int t_rank,
              double *bufSend, double *bufRecv, MPI_Status *status, int chunk )
{
  int myrank;
  int xdim = ir-il+1;
  int ydim = jt-jb+1;
  chunk = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //Send to left &  recieve from right

  if (l_rank != MPI_PROC_NULL)
  {
    for (int j = 1; j <=ydim; ++j)
    {
      bufSend[j-1] = P[1][j];
    }
    MPI_Send(bufSend, ydim, MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD);
  }

  if (r_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, ydim, MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);
    for (int j = 1; j <=ydim; ++j)
    {
      P[xdim+1][j] = bufRecv[j-1];
    }
  }

  // send to right & recieve from left
  if (r_rank != MPI_PROC_NULL)
  {
    for (int j = 1; j <=ydim; ++j)
    {
      bufSend[j-1] = P[xdim][j];
    }
    MPI_Send(bufSend, ydim, MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD);
  }
  if (l_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, ydim, MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
    for (int j = 1; j <=ydim; ++j)
    {
      P[0][j] = bufRecv[j-1];
    }
  }
  //send to top recieve from bottom
  if (t_rank != MPI_PROC_NULL)
  {
    for (int i = 1; i <=xdim; ++i)
    {
      bufSend[i-1] = P[i][ydim];
    }
    MPI_Send(bufSend, xdim, MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD);
  }
  if (b_rank != MPI_PROC_NULL)
  {
    for (int i = 1; i <=xdim; ++i)
    {
      P[i][0] = bufRecv[i-1];
    }
    MPI_Recv(bufRecv, xdim, MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
  }
  

  ///send to bottom recieve from top
  
  if (b_rank != MPI_PROC_NULL)
  {
    for (int i = 1; i <=xdim; ++i)
    {
      bufSend[i-1] = P[i][1];
    }
    MPI_Send(bufSend,xdim, MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD);
  }
  if (t_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, xdim, MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);
    for (int i = 1; i <=xdim; ++i)
    {
      P[i][ydim+1] = bufRecv[i-1];
    }
  }

}

/*
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
  int xdim = ir-il+1;
  int ydim = jt-jb+1;
  chunk = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

 ///////////// Velocity U ////////////
  //Send to left &  recieve from right
  for (int j = 1; j <=ydim; ++j)
  {
    bufSend[j-1] = U[1][j];
  }
  if (l_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, ydim, MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD);
  }
  if (r_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, ydim, MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int j = 1; j <=ydim; ++j)
  {
    U[xdim][j] = bufRecv[j-1];
  }

  // send to right & recieve from left
  for (int j = 1; j <=ydim; ++j)
  {
    bufSend[j-1] = U[xdim-1][j];
  }
  if (r_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, ydim, MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD);
  }
  if (l_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, ydim, MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int j = 1; j <=ydim; ++j)
  {
    U[0][j] = bufRecv[j-1];
  }

  //send to top recieve from bottom
  for (int i = 1; i <=xdim+1; ++i)
  {
    bufSend[i-1] = U[i][ydim];
  }
  if (t_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, (xdim+1), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD);
  }
  if (b_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, (xdim+1), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int i = 1; i <=xdim+1; ++i)
  {
    U[i][0] = bufRecv[i-1];
  }

  ///send to bottom recieve from top
  for (int i = 1; i <=xdim+1; ++i)
  {
    bufSend[i-1] = U[i][1];
  }
  if (b_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, (xdim+1), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD);
      }
  if (t_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, (xdim+1), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int i = 1; i <=xdim+1; ++i)
  {
    U[i][ydim+1] = bufRecv[i-1];
  }


  ///////////// Velocity V///////////////
  //Send to left &  recieve from right
  for (int j = 1; j <=ydim+1; ++j)
  {
    bufSend[j-1] = V[1][j];
  }
  if (l_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, (ydim+1), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD);
  }
  if (r_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, (ydim+1), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int j = 1; j <=ydim+1; ++j)
  {
    V[xdim+1][j] = bufRecv[j-1];
  }

  // send to right & recieve from left
  for (int j = 1; j <=ydim+1; ++j)
  {
    bufSend[j-1] = V[xdim][j];
  }
  if (r_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, (ydim+1), MPI_DOUBLE, r_rank, chunk, MPI_COMM_WORLD);
  }
  if (l_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, (ydim+1), MPI_DOUBLE, l_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int j = 1; j <=ydim+1; ++j)
  {
    V[0][j] = bufRecv[j-1];
  }

  //send to top recieve from bottom
  for (int i = 1; i <=xdim; ++i)
  {
    bufSend[i-1] = V[i][ydim-1];
  }
  if (t_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, (xdim), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD);
  }
  if (b_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, (xdim), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int i = 1; i <=xdim; ++i)
  {
    V[i][0] = bufRecv[i-1];
  }
  ///send to bottom recieve from top
  for (int i = 1; i <=xdim; ++i)
  {
    bufSend[i-1] = V[i][1];
  }
  if (b_rank != MPI_PROC_NULL)
  {
    MPI_Send(bufSend, (xdim), MPI_DOUBLE, b_rank, chunk, MPI_COMM_WORLD);
  }
  if (t_rank != MPI_PROC_NULL)
  {
    MPI_Recv(bufRecv, (xdim), MPI_DOUBLE, t_rank, chunk, MPI_COMM_WORLD, status);
  }
  for (int i = 1; i <=xdim; ++i)
  {
    V[i][ydim] = bufRecv[i-1];
  }

}
*/

void uv_comm (	double**U, double**V,
				int il,int ir,
				int jb, int jt,
				int rank_l, int rank_r,
				int rank_b, int rank_t,
				double *bufSend, double *bufRecv,
				MPI_Status *status, int chunk)
{
	int size;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;

	/* Send to the left, receive from the right */
	/* Send to the right, receive from the left */
	if(rank_l != MPI_PROC_NULL || rank_r != MPI_PROC_NULL)	{
		size = (2 * y_dim) + 1;

		if(rank_l != MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Perform both left-right transfers */
		{

			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy left values to send */
			for(int j = 1; j <= y_dim; j++)
			{
				bufSend[j - 1] = U[1][j];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				bufSend[j - 1] = V[1][j - y_dim];
			}

			/* Send left values, receive right values */
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_l, 1, bufRecv, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);

			/* Copy received right values */
			for(int j = 1; j <= y_dim; j++)
			{
               			 U[x_dim][j] = bufRecv[j - 1];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				V[x_dim + 1][j - y_dim] = bufRecv[j - 1];
			}

			/* Copy right values to send */
			for(int j = 1; j <= y_dim; j++)
			{
				bufSend[j - 1] = U[x_dim-1][j];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				bufSend[j - 1] = V[x_dim][j - y_dim];
			}

			/* Send right values, receive left values */
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_r, 1, bufRecv, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);

			/* Copy received left values */
			for(int j = 1; j <= y_dim; j++)
			{
				U[0][j] = bufRecv[j - 1];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				V[0][j - y_dim] = bufRecv[j - 1];
			}

			free(bufSend);
			free(bufRecv);
		}
		else if(rank_l == MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Receive data from right and  send data to right */
		{

			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			/* Receive right values */
			MPI_Recv(bufRecv, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);

			/* Copy received right values */
			for(int j = 1; j <= y_dim; j++)
			{
				U[x_dim][j] = bufRecv[j - 1];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				V[x_dim + 1][j - y_dim] = bufRecv[j - 1];
			}

			/* Copy right values to send */
			for(int j = 1; j <= y_dim; j++)
			{
				bufSend[j - 1] = U[x_dim-1][j];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				bufSend[j - 1] = V[x_dim][j - y_dim];
			}

			/* Send right values */
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
			free(bufRecv);
			free(bufSend);
		}
		else if(rank_l != MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /*Send data to left and Receive data from left*/
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			/* Copy left values to send */
			for(int j = 1; j <= y_dim; j++)
			{
				bufSend[j - 1] = U[1][j];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				bufSend[j - 1] = V[1][j - y_dim];
			}

			/* Send left values */
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);

			/* Receive left values */
			MPI_Recv(bufRecv, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);

			/* Copy received left values */
			for(int j = 1; j <= y_dim; j++)
			{
				U[0][j] = bufRecv[j - 1];
			}
			for(int j = (y_dim + 1); j <= size; j++)
			{
				V[0][j - y_dim] = bufRecv[j - 1];
			}
			free(bufRecv);
			free(bufSend);
		}
	}

	/* Send to the top, receive from the bottom */
	/* Send to the bottom, receive from the top */
	if(rank_t != MPI_PROC_NULL || rank_b != MPI_PROC_NULL)
	{
		size = (2 * x_dim) + 1;

		if(rank_t != MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /* Perform both top-bottom transfers */
		{

			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy top values to send */
			for(int i = 1; i <= x_dim; i++)
			{
				bufSend[i - 1] = V[i][y_dim-1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				bufSend[i - 1] = U[i - x_dim][y_dim];
			}

			/* Send values to top, receive values from bottom*/
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_t, 1, bufRecv, size, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);

			/* Copy received bottom values */
			for(int i = 1; i <= x_dim; i++)
			{
				V[i][0] = bufRecv[i - 1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				U[i - x_dim][0] = bufRecv[i - 1];
			}

			/* Copy bottom values to send */
			for(int i = 1; i <= x_dim; i++)
			{
				bufSend[i - 1] = V[i][1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				bufSend[i - 1] = U[i - x_dim][1];
			}

			/* Send bottom values, receive top values */
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_b, 1, bufRecv, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);

			/* Copy received top values */
			for(int i = 1; i <= x_dim; i++)
			{
				V[i][y_dim] = bufRecv[i - 1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				U[i - x_dim][y_dim + 1] = bufRecv[i - 1];
			}

			free(bufSend);
			free(bufRecv);
		}
		else if(rank_t == MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /*Receive data from bottom and  send data to bottom*/
		{

			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			/* Receive bottom values */
			MPI_Recv(bufRecv, size, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);

			/* Copy received bottom values */
			for(int i = 1; i <= x_dim; i++)
			{
				V[i][0] = bufRecv[i - 1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				U[i - x_dim][0] = bufRecv[i - 1];
			}

			/* Copy bottom values to send */
			for(int i = 1; i <= x_dim; i++)
			{
				bufSend[i - 1] = V[i][1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				bufSend[i - 1] = U[i - x_dim][1];
			}

			/* Send bottom values */
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);
			free(bufRecv);
			free(bufSend);
		}
		else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /*Send data to top and receive data from top */
		{
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));
			/* Copy top values to send */
			for(int i = 1; i <= x_dim; i++)
			{
				bufSend[i - 1] = V[i][y_dim-1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				bufSend[i - 1] = U[i - x_dim][y_dim];
			}

			/* Send top values */
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);

			/* Receive top values */
			MPI_Recv(bufRecv, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);

			/* Copy received top values */
			for(int i = 1; i <= x_dim; i++)
			{
				V[i][y_dim] = bufRecv[i - 1];
			}
			for(int i = (x_dim + 1); i <= size; i++)
			{
				U[i - x_dim][y_dim + 1] = bufRecv[i - 1];
			}
			free(bufRecv);
			free(bufSend);
		}
	}


}