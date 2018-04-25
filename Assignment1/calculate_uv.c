/*
 * uvp.c
 *
 *  Created on: Apr 20, 2018
 *      Author: abinav
 */

#include <uvp.h>

void calculate_uv(double dt,double dx,double dy,int imax, int jmax, double**U, double**V,double**F,double**G,double **P)
{
	for (int i = 1; i<= imax;++i)
	{
		for (int j=1; j<=jmax;++j)
		{
			//update the U and V components of velocity
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j]-P[i][j]);
			V[i][j] = F[i][j] - (dt/dy)*(P[i][j+1] -P[i][j]);
		}
	}
}
