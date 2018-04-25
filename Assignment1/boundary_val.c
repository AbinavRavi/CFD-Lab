/*
 * boundary_val.c
 *
 *  Created on: Apr 20, 2018
 *      Author: abinav
 */

#include <boundary_val.h>

void boundaryvalues(int imax,int jmax,double **U,double **V)
{
	for (int j =1; j<=jmax;++j)
	{
		//Horizontal velocities on Vertical Boundaries
		U[0,j] = 0;
		U[imax,j] = 0;

		//Vertical Boundary conditions
		 V[0][j]       = - V[1][j];
		 V[imax+1][j]  = - V[imax][j];


	}

	for(int i = 0; i<=imax; ++i)
	{
		//Vertical Velocities on Horizontal Boundaries
		 V[i][0]       = 0;
		 V[i][jmax]    = 0;

		 //Horizontal Boundary Conditions
		 U[i][0]       = - U[i][1];
		 U[i][jmax+1]  = 2.0 - U[i][jmax];

	}
}

