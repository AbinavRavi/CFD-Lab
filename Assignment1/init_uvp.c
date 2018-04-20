#include "helper.h"


void init_uvp(double UI, double VI, double PI, int imax, int jmax, double** U, double** V, double** P){
	init_matrix(U ,1,imax ,1, jmax, UI);
	init_matrix(V ,1,imax ,1, jmax, VI);
	init_matrix(P ,1,imax ,1, jmax, PI);
}
