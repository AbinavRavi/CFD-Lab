#include "helper.h"

void init_uvp(double UI, double VI, double PI, int imax, int jmax, double** U, double** V, double** P){
	init_matrix(U ,0,imax+1 ,0, jmax+1, UI);
	//printf("Debug: init_matrix called inside init_uvp \n");
	init_matrix(V ,0,imax+1 ,0, jmax+1, VI);
	init_matrix(P ,0,imax+1 ,0, jmax+1, PI);
}
