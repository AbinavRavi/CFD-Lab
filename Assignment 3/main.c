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

  int iproc = 2;
  int jproc = 3;
  int imax = 300;
  int jmax = 300;
  int num_proc;
  int *myrank = (int *)malloc(sizeof(int));
  int *il = (int *)malloc(sizeof(int));
  int *ir = (int *)malloc(sizeof(int));
  int *jb = (int *)malloc(sizeof(int));
  int *jt = (int *)malloc(sizeof(int));
  int *rank_l = (int *)malloc(sizeof(int));
  int *rank_r = (int *)malloc(sizeof(int));
  int *rank_b = (int *)malloc(sizeof(int));
  int *rankt = (int *)malloc(sizeof(int));
  int *omg_i = (int *)malloc(sizeof(int));
  int *omg_j = (int *)malloc(sizeof(int));

  MPI_Init(&argn, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, myrank);


  //Initializing the parallel parameters to broadcast
  init_parallel(iproc, jproc, imax, jmax, myrank, il, ir, jb, jt, rank_l,
                rank_r, rank_b, rankt, omg_i, omg_j, num_proc);


  
  printf("This is process %d \n", *myrank);


  MPI_Finalize();
  return 0;
}
