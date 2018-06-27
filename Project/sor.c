#include "sor.h"
#include <math.h>
#include"boundary_val.h"
#include <stdio.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **flag)
{
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

 //set pressure at boundary- no-slip/free-slip/inflow/outflow
for(int i = 0; i<imax; ++i)
{
	for(int j = 0; j<jmax; ++j)
	{
		if ( B_O(flag[i][j]) )  P[i][j] = P[i+1][j];

		if ( B_W(flag[i][j]) )  P[i][j] = P[i-1][j];

		if ( B_N(flag[i][j]) )  P[i][j] = P[i][j+1];

		if ( B_S(flag[i][j]) )  P[i][j] = P[i][j-1];

		if ( B_NO(flag[i][j]) ) P[i][j] = (P[i][j+1] + P[i+1][j])/2;

		if ( B_NW(flag[i][j]) ) P[i][j] = (P[i][j+1] + P[i-1][j])/2;

		if ( B_SO(flag[i][j]) ) P[i][j] = (P[i][j-1] + P[i+1][j])/2;

		if ( B_SW(flag[i][j]) ) P[i][j] = (P[i][j-1] + P[i-1][j])/2;

		if (flag[i][j]&(1<<3) ) P[i][j] = 0;

		if (flag[i][j]&(1<<4) ) P[i][j] = P[i+1][j];
	}
}

  /* SOR iteration */
  for(i = 0; i < imax; i++) {
    for(j = 0; j< jmax; j++) {

	if(flag[i][j]&(1<<0)){

      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*( (P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
	}
    }
  }

  /* compute the residual */
  rloc = 0; int num_fluid_elem=0;

  for(i = 0; i < imax; i++) {
    for(j = 0; j < jmax; j++) {

	if(flag[i][j]&(1<<0)){

      num_fluid_elem++;
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
	}
    }
  }
  rloc = rloc/(num_fluid_elem);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;

//printf("SOR step calculated \n");

}

