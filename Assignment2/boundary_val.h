#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

#include<math.h>
/**
 * The boundary values of the problem are set.
 */
void spec_boundary_val(int imax,int jmax,double **U,double **V,int **flag);

void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int **flag
);


	int B_O(int flag){
	return flag & (1<<8) && ~( flag & ((1<<5) | (1<<6)) );
	}

	int B_W(int flag){
	return flag & (1<<7) && ~( flag & ((1<<5) | (1<<6)) );
	}

	int B_N(int flag){
	return flag & (1<<5) && ~( flag & ((1<<7) | (1<<8)) );
	}

	int B_S(int flag){
	return flag & (1<<6) && ~( flag & ((1<<7) | (1<<8)) );
	}

	int B_NO(int flag){
	return flag & (1<<8) &&  flag & (1<<5);
	}

	int B_NW(int flag){
	return flag & (1<<7) &&  flag & (1<<5);
	}

	int B_SO(int flag){
	return flag & (1<<8) &&  flag & (1<<6); 
	}

	int B_SW(int flag){
	return flag & (1<<7) &&  flag & (1<<6);
	}

#endif
