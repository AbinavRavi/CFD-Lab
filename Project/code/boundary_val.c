#include "boundary_val.h"
#include <stdio.h>

void spec_boundary_val(int imax,int jmax,double **U,double **V,int **flag)
{
	for(int i = 0; i<imax; ++i)
	{
		for(int j = 0; j<jmax; ++j)
		{
		if(flag[i][j]&(1<<4))
			{ U[i][j]=1; V[i][j] = 0; V[i][j-1] = 0;}
		}
	}

}
	
int B_O(int flag){
	return ( (flag & (1<<8)) && ~( flag & ((1<<5) | (1<<6)) ) );
}

int B_W(int flag){
	return ( (flag & (1<<7)) && ~( flag & ((1<<5) | (1<<6)) ) );
}

int B_N(int flag){
	return ( (flag & (1<<5)) && ~( flag & ((1<<7) | (1<<8)) ) );
}

int B_S(int flag){
	return ( (flag & (1<<6)) && ~( flag & ((1<<7) | (1<<8)) ) );
}

int B_NO(int flag){
	return ( (flag & (1<<8)) &&  (flag & (1<<5)) );
}

int B_NW(int flag){
	return ( (flag & (1<<7)) &&  (flag & (1<<5)) );
}

int B_SO(int flag){
	return ( (flag & (1<<8)) &&  (flag & (1<<6)) ); 
}

int B_SW(int flag){
	return ( (flag & (1<<7)) &&  (flag & (1<<6)) );
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

	//spec_boundary_val( imax, jmax, U, V,flag);

	//printf("Boundary values for u,v set. \n");

}


