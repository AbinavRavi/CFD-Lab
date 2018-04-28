#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>


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
    /* for pressure per time step */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;           /* time for output */

    //Read and assign the parameter values from file
    int pseudo;
	pseudo = read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end,
				&xlength, &ylength, &dt, &dx, &dy, &imax, &jmax,
                                &alpha, &omg, &tau, &itermax, &eps, &dt_value);

	pseudo++;
	//printf("Debug: jmax:  %d\n", jmax);
	//printf("Debug: file read \n");
    //Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
    double **P = matrix(0, imax+1, 0, jmax+1);
    double **U = matrix(0, imax+1, 0, jmax+1);
    double **V = matrix(0, imax+1, 0, jmax+1);
    double **F = matrix(0, imax+1, 0, jmax+1);
    double **G = matrix(0, imax+1, 0, jmax+1);
    double **RS = matrix(0, imax+1, 0, jmax+1);
									//printf("Debug: Matrix allocated on heap \n");
    //Initialize the U, V and P
    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
									//printf("Debug: init_uvp called \n");
	int n = 0; //
	double t = 0;
	int n1=0;
//main loop
    while (t < t_end) {
//calling the calculate_dt function for the dt values
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
  // used for testing  printf("%f \n",dt);
  //calling the boundaryvalues function from boundary_val.c
    boundaryvalues(imax,jmax,U,V);
	//calling the function calculate_fg from uvp.c
    calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);
	//calling teh function calculate_rs from uvp.c
    calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);
	//initializing initial iteration to zero
	int it = 0;
  //initializing the residual value with a high number in random
	double res = 10.0;
  //do while loop is ran so as to execute the condition atleast once in case of while condition not being satisfied
    do {
      //calling the sor Function from sor.c
    	sor(omg,dx,dy,imax,jmax,P,RS,&res);

	++it;

    	} while(it<itermax && res>eps);
//calling the calculate_uv function from uvp.c
  	calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);
  	if (t >= n1*dt_value)
  	{
      //writing the output values
   		write_vtkFile("solution", n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
		printf("%f SECONDS COMPLETED \n",n1*dt_value);
    		n1=n1+ 1;
    		continue;
  	}
    t =t+ dt; //increase in timestep dt
    n = n+ 1; 
    }

    //Free memory
    free_matrix( P, 0, imax, 0, jmax);
    free_matrix( U, 0, imax, 0, jmax);
    free_matrix( V, 0, imax, 0, jmax);
    free_matrix( F, 0, imax, 0, jmax);
    free_matrix( G, 0, imax, 0, jmax);
    free_matrix(RS, 0, imax, 0, jmax);

  return -1;

}
