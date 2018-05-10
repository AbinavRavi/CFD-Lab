#include "boundary_val.h"
#include <stdio.h>

void spec_boundary_val(int imax,int jmax,double **U,double **V,int **flag){
	for(int i = 0; i<imax; ++i)
	{
		for(int j = 0; j<jmax; ++j)
		{
		if(flag[i][j]&(1<<4))
			{ U[i][j]=1; V[i][j] = 0; V[i][j-1] = 0;}
		}
	}

}

void boundaryvalues(int imax,int jmax,double **U,double **V,int **flag)
{
	for(int i = 0; i<imax; ++i)
	{
	for(int j = 0; j<jmax; ++j)
	{

	switch(flag[i][j]& ((1<<1)|(1<<2)|(1<<3)) )
{
	case 1<<1://No slip conditions
	if ( B_O(flag[i][j]) )
	{
	U[i][j] = 0;
	V[i][j-1] = -V[i+1][j-1];	
	V[i][j] = -V[i+1][j];
	}

	if ( B_W(flag[i][j]) ) 
	{
	U[i-1][j] = 0;
	V[i][j-1] = -V[i-1][j-1];	
	V[i][j] = -V[i-1][j];
	}

	if ( B_N(flag[i][j]) ) 
	{
	V[i][j] = 0;
	U[i-1][j] = -U[i-1][j+1];
	U[i][j] = -U[i][j+1];	
	}

	if ( B_S(flag[i][j]) )
	{
	V[i][j-1] = 0;
	U[i-1][j] = -U[i-1][j-1];
	U[i][j] = -U[i][j-1];	
	}

	if ( B_NO(flag[i][j])  )
	{
	U[i][j] = 0;
	U[i-1][j] = -U[i-1][j+1];
	V[i][j] = 0;
	V[i][j-1] = -V[i+1][j-1];
	}

	if ( B_NW(flag[i][j])  )
	{
	U[i-1][j] = 0;
	U[i][j] = - U[i][j+1];
	V[i][j] = 0;
	V[i][j-1] = -V[i-1][j-1];
	} 

	if ( B_SO(flag[i][j]) ){
	U[i][j]=0;
	U[i-1][j] = -U[i-1][j-1];
	V[i][j-1]=0;
	V[i][j] = -V[i+1][j];
	}

	if ( B_SW(flag[i][j]) ){
	U[i-1][j] = 0;
	U[i][j] = -U[i][j-1];
	V[i][j-1] = 0;
	V[i][j] = -V[i-1][j];
	}
	break;

	case 1<<2://Free Slip conditions
	if ( B_O(flag[i][j]) )
	{
	U[i][j] = 0;	
	V[i][j] = V[i+1][j];
	V[i][j-1] = V[i+1][j-1];
	}

	if ( B_W(flag[i][j]) ) 
	{
	U[i-1][j] = 0;
	V[i][j-1] = V[i-1][j-1];	
	V[i][j] = V[i-1][j];
	}

	if ( B_N(flag[i][j]) ) 
	{
	V[i][j] = 0;
	U[i-1][j] = U[i-1][j+1];
	U[i][j] = U[i][j+1];	
	}

	if ( B_S(flag[i][j]) )
	{
	V[i][j-1] = 0;
	U[i-1][j] = U[i-1][j-1];
	U[i][j] = U[i][j-1];	
	}

	if ( B_NO(flag[i][j])  )
	{
	U[i][j] = 0;
	U[i-1][j] = U[i-1][j+1];
	V[i][j] = 0;
	V[i][j-1] = V[i+1][j-1];
	}

	if ( B_NW(flag[i][j])  )
	{
	U[i-1][j] = 0;
	U[i][j] = U[i][j+1];
	V[i][j] = 0;
	V[i][j-1] = V[i-1][j-1];
	} 

	if ( B_SO(flag[i][j]) ){
	U[i][j]=0;
	U[i-1][j] = U[i-1][j-1];
	V[i][j-1]=0;
	V[i][j] = V[i+1][j];
	}

	if ( B_SW(flag[i][j]) ){
	U[i-1][j] = 0;
	U[i][j] = U[i][j-1];
	V[i][j-1] = 0;
	V[i][j] = V[i-1][j];
	}	
	break;

	case 1<<3://Assuming outflow happens in x direction.
	U[i][j] = U[i-1][j];
	V[i][j] = V[i-1][j];
	V[i][j-1] = V[i-1][j-1];
	break;
}
	}
	}

	spec_boundary_val( imax, jmax, U, V,flag);

	printf("Boundary values for u,v set. \n");

}


