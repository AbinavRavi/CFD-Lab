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

    printf("Assignment-2, Group D \n");
    printf("Please select the problem from the list below by typing 1-5 \n");
    printf("1. Karman Vortex Street \n");
    printf("2. Flow over a Step \n");
    printf("3. Natural Convection \n");
    printf("4. Fluid Trap \n");
    printf("5. rayleigh-Benard Convection \n");
    int select;
    scanf("%d",&select);
    //location of input file
    const char* filename = "0";
    switch(select)
    {
    case 1:
    filename = "karman_vortex.dat";
    break;
    case 2:
    filename = "step_flow.dat";
    break;
    case 3:
    filename = "natural_convection.dat";
    break;
    case 4:
    filename = "fluid_trap.dat";
    break;
    case 5:
    filename = "RB_convection.dat";
    break;
    }

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
    double Pr;
    double TI;
    double T_h;
    double T_c;
    double beta;
    char *problem = "step_flow";
    char *geometry = "step_flow.pgm";
    //Read and assign the parameter values from file
    read_parameters(filename, &imax, &jmax, &xlength, &ylength, 
			&dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax,
			&GX, &GY, &Re, &Pr, &UI, &VI, &PI, &TI, &T_h, &T_c, &beta, &dx, &dy);

    
    int include_temp = 1;
    if(((select==1)||(select==2)))
	{
		if( (Pr!=0)||(TI!=0)||(T_h!=0)||(T_c!=0)||(beta!=0) ){
		char szBuff[80];
        sprintf( szBuff, "Input file incompatible. Please check .dat file. \n");
        ERROR( szBuff );
		}
		else  include_temp = 0;
	}

    //Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
    printf("PROGRESS: Starting matrix allocation... \n");
    double **P = matrix(0, imax-1, 0, jmax-1);
    double **U = matrix(0, imax-1, 0, jmax-1);
    double **V = matrix(0, imax-1, 0, jmax-1);
    double **F = matrix(0, imax-1, 0, jmax-1);
    double **G = matrix(0, imax-1, 0, jmax-1);
    double **RS = matrix(0, imax-1, 0, jmax-1);
    int **flag = imatrix(0, imax-1, 0, jmax-1);
	if(include_temp) double **temp = matrix(0, imax-1, 0, jmax-1);
    printf("PROGRESS: Matrices allocated on heap... \n \n");

    //Initilize flags
    init_flag(problem,geometry, imax, jmax, flag);

    //Initialize the U, V and P
    if(include_temp)
	{
		init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, temp, flag);
	}
	else
	{
		init_uvp(UI, VI, PI, imax, jmax, U, V, P, flag);
	}


    printf("PROGRESS: Starting the flow simulation...\n");
    double t=0; int n=0; int n1=0;
    while (t < t_end) {
        const char* is_converged = "Yes";
	calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
   	printf("t = %f ,dt = %f, ",t,dt);
							
    	boundaryvalues(imax, jmax, U, V, flag);

    	spec_boundary_val(imax, jmax, U, V, flag);
			
    	calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G,flag);
													
    	calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);
													
	int it = 0;
	double res = 10.0;

    	do {

    		sor(omg,dx,dy,imax,jmax,P,RS,&res,flag);
    	
		++it;

    	} while(it<itermax && res>eps);
	printf("SOR itertions = %d ,residual = %f \n", it-1, res);
	if((it==itermax)&&(res>eps)){
		printf("WARNING: Iteration limit reached before convergence. \n");
		is_converged = "No";
	}
  	write_sim_log(problem, t, dt, n, it-1, res, is_converged);

	calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);
	
  	if (t >= n1*dt_value)
  	{
   		write_vtkFile("solution", n,xlength,ylength,imax-2,jmax-2,dx,dy,U,V,P);
		printf("writing result at %f SECONDS COMPLETED \n",n1*dt_value);
    		n1=n1+ 1;
    		continue;
  	}
    	t =t+ dt;
    	n = n+ 1;
    }
    printf("PROGRESS: flow simulation completed...\n");

    printf("PROGRESS: Freeing allocated memory...\n");
    //Free memory
    free_matrix( P, 0, imax, 0, jmax);
    free_matrix( U, 0, imax, 0, jmax);
    free_matrix( V, 0, imax, 0, jmax);
    free_matrix( F, 0, imax, 0, jmax);
    free_matrix( G, 0, imax, 0, jmax);
    free_matrix(RS, 0, imax, 0, jmax);
    free_imatrix(flag, 0, imax, 0, jmax);
    printf("PROGRESS: allocated memory released...\n");
  return -1;
    
}


/*Things to do:
read parameter
write temp and obstructions to vtk filename
assert function
*/
