#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>
#include "sor.h"
#include "uvp.h"
#include "boundary_val.h"


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

  int n = 0; //
  int t = 0;
  int n1=1;
  while (t < tend) {
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
    boundaryvalues(imax,jmax,U,V);
    calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);
    calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);
	int it = 0;
	double res = 0.0;
  do {
    sor(omg,dx,dy,imax,jmax,P,RS,&res);
    ++it;
  } while(it<itmax && res>eps);
  calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);
  if (t >= n1*dt_end)
  {
    write_vtkFile("solution", n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
    n1+=1;
    break;
  }
  t +=dt;
  n = n+1;
  }
  return -1;
}
