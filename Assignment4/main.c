#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "/home/snowcrash/CSESS18/CFDLab/precice-1.1.1/src/precice/adapters/c/SolverInterfaceC.h"
#include "/home/snowcrash/CSESS18/CFDLab/precice-1.1.1/src/precice/adapters/c/Constants.h"
#include"precice_adapter.h"


int main(int argn, char** args){

	printf("Start of Run... \n");
			printf("Assignment-4, Group D \n");
			printf("Please select the problem from the list below by typing 1-4 \n");
			printf("P1. Forced Convection over a heated plate \n");
			printf("P2. Natural Convection with heat conducting walls \n");
			printf("P3. 2D Heat Exchanger - Fluid1 \n");
			printf("P4. 2D Heat Exchanger - Fluid2 \n");
			int select;
			char* geometry = (char*)(malloc(sizeof(char)*200));
			char* problem = (char*)(malloc(sizeof(char)*2));
			scanf("%d",&select);
			//select problem
			const char* filename = "0";
			switch(select)
			{
			case 1:
			filename = "configs/heated-plate.dat";
			break;

			case 2:
			filename = "configs/convection.dat";
			break;

			case 3:
			filename = "configs/F1-heat-exchange.dat";
			break;

			case 4:
			filename = "configs/F2-heat-exchange.dat";
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
	double x_origin;
	double y_origin;
    double dt ;                /* time step */
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
    double Pr;
    double TI;
    double T_h;
    double T_c;
    double beta;

	char *precice_config = (char*)(malloc(sizeof(char)*200));
	char *participant_name = (char*)(malloc(sizeof(char)*200));
	char *mesh_name = (char*)(malloc(sizeof(char)*200));
	char *read_data_name = (char*)(malloc(sizeof(char)*200));
	char *write_data_name = (char*)(malloc(sizeof(char)*200));

	const char* coric = precicec_actionReadIterationCheckpoint();
	const char* cowic = precicec_actionWriteIterationCheckpoint();

    //Read and assign the parameter values from file
    read_parameters(filename, &imax, &jmax, &xlength, &ylength,&x_origin, &y_origin,
			&dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax,
			&GX, &GY, &Re, &Pr, &UI, &VI, &PI, &TI, &beta, &dx, &dy, problem, geometry,
			precice_config, participant_name, mesh_name, read_data_name, write_data_name);
    
    //Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
    printf("PROGRESS: Starting matrix allocation... \n");
    double **P = matrix(0, imax-1, 0, jmax-1);
    double **U = matrix(0, imax-1, 0, jmax-1);
    double **V = matrix(0, imax-1, 0, jmax-1);
    double **U_cp = matrix(0, imax-1, 0, jmax-1);
    double **V_cp = matrix(0, imax-1, 0, jmax-1);
    double **F = matrix(0, imax-1, 0, jmax-1);
    double **G = matrix(0, imax-1, 0, jmax-1);
    double **RS = matrix(0, imax-1, 0, jmax-1);
    int **flag = imatrix(0, imax-1, 0, jmax-1);
    double **T = matrix(0, imax-1, 0, jmax-1);
	double **T_cp = matrix(0, imax-1, 0, jmax-1);
    double **T1 = matrix(0, imax-1, 0, jmax-1);

    printf("PROGRESS: Matrices allocated on heap... \n \n");

	//Initialize precice
	precicec_createSolverInterface(participant_name, precice_config, 0, 1);
	int dim = precicec_getDimensions();
	
	// define coupling mesh
	int meshID = precicec_getMeshID(mesh_name);
	int num_coupling_cells = 0;

    //Initilize flags
    init_flag(problem,geometry, imax, jmax, flag, &num_coupling_cells);

    //Initialize the U, V and P
	init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag);

	// get coupling cell ids
	int* vertexIDs = precice_set_interface_vertices(imax, jmax, dx, dy, x_origin, y_origin,
                                					num_coupling_cells, meshID, flag); 

	// define Dirichlet part of coupling written by this solver
	int temperatureID = precicec_getDataID(write_data_name, meshID);

	double* temperatureCoupled = (double*) malloc(sizeof(double)*num_coupling_cells);

	// define Neumann part of coupling read by this solver
	int heatFluxID = precicec_getDataID(read_data_name, meshID);
	double* heatFluxCoupled = (double*) malloc(sizeof(double) * num_coupling_cells);

	// call precicec_initialize()
	double precice_dt = precicec_initialize();

	// initialize data at coupling interface
	precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled,
							vertexIDs, temperatureID, T, flag);
	// synchronize with OpenFOAM
	precicec_initialize_data();
	// read heatfluxCoupled
	precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatFluxCoupled); 

	//Make solution folder
	struct stat st = {0};
	char sol_folder[80];
	sprintf( sol_folder,"Fluid_%s",problem);
	if (stat(sol_folder, &st) == -1) {
    		mkdir(sol_folder, 0700);
	}
	char sol_directory[80];
	sprintf( sol_directory,"Fluid_%s/sol", problem);
	//create log file
	char LogFileName[80];
 	FILE *fp_log = NULL;
	sprintf( LogFileName, "%s.log", problem );
	fp_log = fopen( LogFileName, "w");
	fprintf(fp_log, "It.no.|   Time    |time step |SOR iterations | residual | SOR converged \n");


    printf("PROGRESS: Starting the flow simulation...\n");
    double t=0; int n=0; int n1=0;
	double time_cp;

	while (precicec_isCouplingOngoing() ) {

		if(precicec_isActionRequired(cowic)){
			write_checkpoint(t, U, V, T, &time_cp, U_cp, V_cp, T_cp, imax, jmax);
			precicec_fulfilledAction(cowic);
		}

        char* is_converged = "Yes";
		//Calculate time step using min of precice_dt and dt
		calculate_dt(Re, tau,&dt,dx,dy,imax,jmax, U, V, Pr);

		dt = fmin(dt, precice_dt);
		printf("t = %f ,dt = %f, ",t,dt);
		//set boundary values for U and V in the fluid domain
		boundaryvalues(imax, jmax, U, V, flag,UI);

		//Set coupling neumann boundary using precice at the coupling interface.
		//this value comes from the adaptor which gives the heat-flux from solid domain computation
		set_coupling_boundary(imax, jmax, dx, dy, heatFluxCoupled, T, flag);

		//calculate the temperature in the fluid
		calculate_temp(T, T1, Pr, Re, imax, jmax, dx, dy, dt, alpha, U, V, flag, TI, UI);

		//set inlet boundary conditions for case-specific problem
    	spec_boundary_val(imax, jmax, U, V, flag,select);
		//calculate F and G

    	calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G,flag, beta, T,UI);
		//Calculate RS

    	calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);

		int it = 0;
		double res = 10.0;
		// SOR iteration to calculate the pressure in fluid domain
    	do {

    		sor(omg,dx,dy,imax,jmax,P,RS,&res,flag,UI);
			++it;

    	} while(it<itermax && res>eps);
		printf("SOR itertions = %d ,residual = %f \n", it-1, res);
		if((it==itermax)&&(res>eps)){
			printf("WARNING: Iteration limit reached before convergence. \n");
			is_converged = "No";
		}
  		fprintf(fp_log, "    %d |  %f | %f |      %d      | %f | %s \n", n, t, dt, it-1, res, is_converged);
		// calculate the velocities U and V
		calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag, UI);

		// write temperatures at the coupling interface to be used by adaptor. This value is given to openfoam
		// to be used for solid domain computation
		precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled,
							vertexIDs, temperatureID, T, flag);
		// get time step from the coupling interface for the solid domain
		precice_dt = precicec_advance(dt);
		// read the heat flux from the solid domain into an array
		precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatFluxCoupled);
		//write data into .vtk files for visualization
		if ((t >= n1*dt_value)&&(t!=0.0))
  		{
			nullify_obstacles(U, V, P, T, flag, imax, jmax);//nullify values at obstacle cells
   			write_vtkFile(sol_directory ,n ,xlength ,ylength, x_origin,y_origin,imax-2 ,jmax-2 ,
							dx ,dy ,U ,V ,P,T);

			printf("writing result at %f seconds \n",n1*dt_value);
    		n1=n1+ 1;
    		continue;
  		}

		if(precicec_isActionRequired(coric)){ // timestep not converged
			restore_checkpoint(&t, U, V, T, time_cp, U_cp, V_cp, T_cp, imax, jmax);
			precicec_fulfilledAction(coric);
		}
		else{ // timestep converged
			t =t+ dt;
			n = n+ 1;
		}

    }
	
	precicec_finalize();

	fclose(fp_log);
    printf("PROGRESS: flow simulation completed...\n \n");

    printf("PROGRESS: Freeing allocated memory...\n");
    //Free memory
    free_matrix( P, 0, imax-1, 0, jmax-1);
    free_matrix( U, 0, imax-1, 0, jmax-1);
    free_matrix( V, 0, imax-1, 0, jmax-1);
    free_matrix( U_cp, 0, imax-1, 0, jmax-1);
    free_matrix( V_cp, 0, imax-1, 0, jmax-1);
    free_matrix( F, 0, imax-1, 0, jmax-1);
    free_matrix( G, 0, imax-1, 0, jmax-1);
    free_matrix(RS, 0, imax-1, 0, jmax-1);
    free_imatrix(flag, 0, imax-1, 0, jmax-1);
	free_matrix(T, 0, imax-1, 0, jmax-1);
	free_matrix(T_cp, 0, imax-1, 0, jmax-1);
	free_matrix(T1, 0, imax-1, 0, jmax-1);
	free(geometry);
	free(problem);
	free(precice_config);
	free(participant_name);
	free(mesh_name);
	free(read_data_name);
	free(write_data_name);
	free(vertexIDs);
	free(temperatureCoupled);
	free(heatFluxCoupled);
    printf("PROGRESS: allocated memory released...\n \n");

	printf("PROGRESS: End of Run.\n");
  return -1;

}
