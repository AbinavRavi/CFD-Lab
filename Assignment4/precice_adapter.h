//
// Created by Andreas Reiser on 27.02.18.
//

#ifndef MT_SOURCE_PRECICE_ADAPTER_H
#define MT_SOURCE_PRECICE_ADAPTER_H

#include<stdlib.h>
#include<stdio.h>
#include "/home/snowcrash/CSESS18/CFDLab/precice-1.1.1/src/precice/adapters/c/SolverInterfaceC.h"


void write_checkpoint(double time, double **U, double **V, double **TEMP, double *time_cp, double **U_cp, double **V_cp,
                      double **TEMP_cp, int imax, int jmax);
void restore_checkpoint(double *time, double **U, double **V, double **TEMP, double time_cp, double **U_cp,
                        double **V_cp,
                        double **TEMP_cp, int imax, int jmax);


int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells,int meshID, int **flag);
void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs,
                               int temperatureID, double **TEMP, int **flag);
void set_coupling_boundary(int imax, int jmax, double dx, double dy, double *heatflux, double **TEMP, int **flag);


#endif //MT_SOURCE_PRECICE_ADAPTER_H
