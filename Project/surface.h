#ifndef SURFACE_H
#define SURFACE_H

#include "helper.h"

struct particle{
    double x, y;
    struct particle *next;
};

struct particleline{
    int length;
    struct particle *Particles;
};

//Surface type 1
int S_O(int flag);

int S_W(int flag);

int S_N(int flag);

int S_S(int flag);

// Surface type 2
int S_NO(int flag);

int S_OS(int flag);

int S_SW(int flag);

int S_WN(int flag);

//Surface type 3
int S_OW(int flag);

int S_NS(int flag);

//Surface type 4
int S_NOS(int flag);

int S_OSW(int flag);

int S_SWN(int flag);

int S_WNO(int flag);

//Surface type 5
int S_NOSW(int flag);

struct particleline *INIT_PARTICLES (int *N, int imax, int jmax, double delx, double dely, int ppc, int **flag, int select);

void insert_particles(struct particle *start, int i, int j, double delx, double dely, int len, int cell_count);

void MARK_CELLS(int **flag, int imax, int jmax, double delx, double dely, int N, struct particleline *Partlines);

void SET_UVP_SURFACE(double **U,double **V, double **P, int **flag, int imax, int jmax, double Re, double delx, double dely, double delt,double GX, double GY);

double U_interp(double **U, double dx, double dy, double x, double y);

double V_interp(double **V, double dx, double dy, double x, double y);

void ADVANCE_PARTICLES(double **U, double **V, double delx, double dely, double delt, int N, struct particleline *Partlines, int **flag, int imax, int jmax);

void FREE_PARTICLELINES(struct particleline *Partlines, int N, int imax, int jmax, double delx, double dely, int ppc, int **flag );

void free_particles(struct particle *start);


#endif