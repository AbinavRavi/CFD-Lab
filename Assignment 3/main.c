#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>
#include"parallel.h"

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

  int iproc;
  int jproc;
  int num_proc;
  int myrank
  int il;
  int ir;
  int jb;
  int jt;
  int l_rank;
  int r_rank;
  int b_rank;
  int t_rank;
  int omg_i;
  int omg_j;

//Read and assign the parameter values from file
	  read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, 
				&xlength, &ylength, &dt, &dx, &dy, &imax, &jmax,
                                &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
  
  MPI_Status *status;

  MPI_Init(&argn, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, myrank);


  //Initializing the parallel parameters to broadcast
  init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &l_rank,
                &r_rank, &b_rank, &t_rank, &omg_i, &omg_j, num_proc);

  
//Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
    double **P = matrix(0, ir-il+1, 0, jt-jb+1);
    double **U = matrix(0, ir-il+1, 0, jt-jb+1);
    double **V = matrix(0, ir-il+1, 0, jt-jb+1);
    double **F = matrix(0, ir-il+1, 0, jt-jb+1);
    double **G = matrix(0, ir-il+1, 0, jt-jb+1);
    double **RS = matrix(0, ir-il+1, 0, jt-jb+1);

//Initialize the U, V and P
    //init_uvp(UI, VI, PI, imax, jmax, U, V, P);
  int n = 0; //
	double t = 0;
	int n1=0;
	int it;
	double res;
  while (t < t_end)
  {
    
    boundaryvalues(imax,jmax,U,V);
													
    calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);
													
    calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);
													
	  it = 0;
	  res = 10.0;

    do {

    	sor(omg,dx,dy,imax,jmax,P,RS,&res);
    	pressure_comm();
//synchronize
	    ++it;
    	} while(it<itermax && res>eps);

  	calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);
  	if (t >= n1*dt_value)
  	{
   		write_vtkFile("solution", n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
		  printf("%f SECONDS COMPLETED \n",n1*dt_value);
    		n1=n1+ 1;
    		continue;
  	}
	  calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
    t =t+ dt;
    n = n+ 1;
    }

    //Free memory
    free_matrix( P, 0, ir-il+1, 0, jt-jb+1);
    free_matrix( U, 0, ir-il+1, 0, jt-jb+1;
    free_matrix( V, 0, ir-il+1, 0, jt-jb+1);
    free_matrix( F, 0, ir-il+1, 0, jt-jb+1);
    free_matrix( G, 0, ir-il+1, 0, jt-jb+1);
    free_matrix(RS, 0, ir-il+1, 0, jt-jb+1);

  Program_Stop();
  return 0;
}
