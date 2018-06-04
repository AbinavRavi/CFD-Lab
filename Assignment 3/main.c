#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include "parallel.h"
#include "mpi.h"
#include <time.h>
/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){

      //location of input file
    const char* filename = "cavity100.dat";

    //define parameter variables
    double Re;                /* reynolds number   */
    double UI;                /* velocity x-direction */
    double VI;                /* velocity y-direction */
    double PI;                /* pressure */
    double GX;                /* gravitation x-direction */
    double GY;                /* gravitation y-direction */
    double t_end;             /* end time */
    double xlength;           /* length of the domain x-dir.*/
    double ylength;           /* length of the domain y-dir.*/
    double dt;                /* time step */
    double dx;                /* length of a cell x-dir. */
    double dy;                /* length of a cell y-dir. */
    int  imax;                /* number of cells x-direction*/
    int  jmax;                /* number of cells y-direction*/
    double alpha;             /* uppwind differencing factor*/
    double omg;               /* relaxation factor */
    double tau;               /* safety factor for time step*/
    int  itermax;             /* max. number of iterations  */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;           /* time for output */

  int iproc = 0;
  int jproc = 0;
  int num_proc = 0;
  int myrank = 0;
  int il = 0;
  int ir = 0;
  int jb = 0;
  int jt = 0;
  int l_rank = 0;
  int r_rank = 0;
  int b_rank = 0;
  int t_rank = 0;
  int omg_i = 0;
  int omg_j = 0;



  MPI_Init(&argn, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//Read and assign the parameter values from file
  int pseudo = 0;
	pseudo = read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end,
				                  &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax,
                          &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);

  pseudo++;
  Programm_Sync("Input file read... \n");
  MPI_Status status;

  //Initializing the parallel parameters to broadcast
  init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &l_rank,
                &r_rank, &b_rank, &t_rank, &omg_i, &omg_j, num_proc);


  Programm_Sync("Parallel parameters initialized... \n");

  int xdim = ir-il+1;
  int ydim = jt-jb+1;

//Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
    double **P = matrix(0, xdim+1, 0, ydim+1);
    double **U = matrix(0, xdim+2, 0, ydim+1);
    double **V = matrix(0, xdim+1, 0, ydim+2);
    double **F = matrix(0, xdim+2, 0, ydim+1);
    double **G = matrix(0, xdim+1, 0, ydim+2);
    double **RS = matrix(0, xdim+1, 0, ydim+1);

//Initialize the U, V and P
  Programm_Sync("U,V, P initializing... \n");
  init_uvp(UI, VI, PI, xdim, ydim, U, V, P);

  int n = 0; //
	double t = 0;
	int n1=0;
	int it;
	double res;
  int buffsize = (xdim)>(ydim)?(xdim+1):(ydim+1);
  double *bufSend = (double*)malloc(sizeof(double)*2*buffsize);
  double *bufRecv = (double*)malloc(sizeof(double)*2*buffsize);
  int chunk = 0;
  Programm_Sync("Starting the simulation... \n");
  
  while (t < t_end)
  {
    boundaryvalues(xdim, ydim, U, V,l_rank, r_rank, b_rank, t_rank );

    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, xdim, ydim, U, V, F, G, l_rank,r_rank, b_rank, t_rank);

    calculate_rs(dt, dx, dy, xdim, ydim, F, G, RS);

	  it = 0;
	  res = 1.0;

    do {
     sor(omg, dx, dy, P, RS, &res, il, ir, jb, jt, l_rank, r_rank, b_rank, t_rank, bufSend, bufRecv, imax,jmax,status);

     MPI_Barrier(MPI_COMM_WORLD);
	    ++it;
    } while(it<itermax&& res>eps);

  	calculate_uv(dt, dx, dy, xdim, ydim, U, V, F, G, P);
    boundaryvalues(xdim, ydim, U, V,l_rank, r_rank, b_rank, t_rank );
    uv_comm(U, V, il, ir, jb, jt, l_rank, r_rank, b_rank, t_rank, bufSend, bufRecv, &status, chunk);

    calculate_dt(Re, tau, &dt, dx, dy, xdim, ydim, U, V);

    if (t >= n1*dt_value)
  	{
		  output_uvp(U, V, P, il, ir, jb, jt, omg_i, omg_j, "Solution", n1);
      printf("%f SECONDS COMPLETED \n",n1*dt_value);
    		n1=n1+ 1;
  	}
    if(myrank ==0){
      printf("t = %f ,dt = %f, Res = %f,iterations=%d \n",t,dt,res,it-1);
    }
    t =t+ dt;
    n = n+ 1;
  }
    Programm_Sync("End of simulation... \n");
    //Free memory
    free_matrix( P, 0, xdim+1, 0, ydim+1);
    free_matrix( U, 0, xdim+2, 0, ydim+1);
    free_matrix( V, 0, xdim+1, 0, ydim+2);
    free_matrix( F, 0, xdim+2, 0, ydim+1);
    free_matrix( G, 0, xdim+1, 0, ydim+2);
    free_matrix(RS, 0, xdim+1, 0, ydim+1);
    free(bufSend);
    free(bufRecv);

    Programm_Stop("Exit");
  return 0;
}
