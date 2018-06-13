#include "precice_adapter.h"


// Check for  coupling cell
int isCouplingCell(int flag)
{
    if (flag&&1<<9!=0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


int *precice_set_interface_vertices(int imax, int jmax, int dx, int dy, double x_origin, double y_origin,
                                    int num_coupling_cells,int meshID, int **flag)
{
    double *vertices = (double*)malloc(sizeof(double)*num_coupling_cells*3);
    int *vertexID = (int*)malloc(sizeof(int)*num_coupling_cells);
    int k = 0;
    //moving over the walls
    //bottom wall
    for(int i = 0; i<imax; i++)
        if (isCouplingCell(flag[i][0])==1) {
            vertices[3*k]= x_origin + (i-0.5)*dx;
            vertices[3*k+1] = y_origin;
            vertices[3*k+2] = 0;
            k++;
        }
    //top wall
    for(int i = 0; i<imax; i++)
        if (isCouplingCell(flag[i][jmax-1])==1) {
            vertices[3*k]= x_origin + (i-0.5)*dx;
            vertices[3*k+1] = y_origin + (jmax-1)*dy;
            vertices[3*k+2] = 0;
            k++;
        }
    //left wall
    for(int j = 0 ; j< jmax; j++)
        if (isCouplingCell(flag[0][j])==1) {
            vertices[3*k]= x_origin;
            vertices[3*k+1] = y_origin + (j-0.5)*dy;
            vertices[3*k+2] = 0;
            k++;
        }
    //right wall
    for(int j = 0 ; j< jmax; j++)
        if (isCouplingCell(flag[imax-1][j])==1) {
            vertices[3*k]= x_origin + (imax-1)*dx;
            vertices[3*k+1] = y_origin + (j-0.5)*dy;
            vertices[3*k+2] = 0;
            k++;
        }

    //inside the domain
    for(int j = 1 ; j< jmax-1; j++)
        for(int i = 0; i<imax; i++)
        {
            //if bottom is fluid
            if ( (isCouplingCell(flag[i][j])==1)&&(flag[i][j-1]&1<<0) ) {
                vertices[3*k]= x_origin + (i-0.5)*dx;
                vertices[3*k+1] = y_origin + (j-1)*dy;
                vertices[3*k+2] = 0;
                k++;
            }
            //if top is fluid
            if ( (isCouplingCell(flag[i][j])==1)&&(flag[i][j+1]&1<<0) ) {
                vertices[3*k]= x_origin + (i-0.5)*dx;
                vertices[3*k+1] = y_origin + (j+1)*dy;
                vertices[3*k+2] = 0;
                k++;
            }
        }
    //assigning values to vertex array
    precicec_setMeshVertices(meshID,num_coupling_cells, vertices ,vertexID);

    return(*vertexID);
    free(vertices);
}


void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs,
                                int temperatureID, double **TEMP, int **flag)
{
    int k=0;    
    //moving over the walls
    //bottom wall
    for(int i = 0; i<imax; i++)
        if (isCouplingCell(flag[i][0])==1) {
            temperature[k] = TEMP[i][1];
            k++;
        }
    //top wall
    for(int i = 0; i<imax; i++)
        if (isCouplingCell(flag[i][jmax-1])==1) {
            temperature[k] = TEMP[i][jmax-2];
            k++;
        }
    //left wall
    for(int j = 0 ; j< jmax; j++)
        if (isCouplingCell(flag[0][j])==1) {
            temperature[k] = TEMP[1][j];
            k++;
        }
    //right wall
    for(int j = 0 ; j< jmax; j++)
        if (isCouplingCell(flag[imax-1][j])==1) {
            temperature[k] = TEMP[imax-2][j];
            k++;
        }

    //inside the domain
    for(int j = 1 ; j< jmax-1; j++)
        for(int i = 0; i<imax; i++)
        {
            //if bottom is fluid
            if ( (isCouplingCell(flag[i][j])==1)&&(flag[i][j-1]&1<<0) ) {
                temperature[k] = TEMP[i][j-1];
                k++;
            }
            //if top is fluid
            if ( (isCouplingCell(flag[i][j])==1)&&(flag[i][j+1]&1<<0) ) {
                temperature[k] = TEMP[i][j+1];
                k++;
            }
        }
    precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
}


void set_coupling_boundary(int imax, int jmax, double dx, double dy, 
                            double *heatflux, double **TEMP, int **flag)
{
    int k=0;    
    //moving over the walls
    //bottom wall
    for(int i = 0; i<imax; i++)
        if (isCouplingCell(flag[i][0])==1) {
            TEMP[i][0] = TEMP[i][1] + dy*heatflux[k];
            k++;
        }
    //top wall
    for(int i = 0; i<imax; i++)
        if (isCouplingCell(flag[i][jmax-1])==1) {
            TEMP[i][jmax-1] = TEMP[i][jmax-2] + dy*heatflux[k];
            k++;
        }
    //left wall
    for(int j = 0 ; j< jmax; j++)
        if (isCouplingCell(flag[0][j])==1) {
            TEMP[0][j] = TEMP[1][j] + dx*heatflux[k];
            k++;
        }
    //right wall
    for(int j = 0 ; j< jmax; j++)
        if (isCouplingCell(flag[imax-1][j])==1) {
            TEMP[imax-1][j] = TEMP[imax-2][j] + dx*heatflux[k];
            k++;
        }

    //inside the domain
    for(int j = 1 ; j< jmax-1; j++)
        for(int i = 0; i<imax; i++)
        {
            //if bottom is fluid
            if ( (isCouplingCell(flag[i][j])==1)&&(flag[i][j-1]&1<<0) ) {
                TEMP[i][j] = TEMP[i][j-1] + dy*heatflux[k];
                k++;
            }
            //if top is fluid
            if ( (isCouplingCell(flag[i][j])==1)&&(flag[i][j+1]&1<<0) ) {
                TEMP[i][j] = TEMP[i][j+1] + dy*heatflux[k];
                k++;
            }
        }
}

/*
void write_checkpoint(double time, double **U, double **V, double **TEMP, double *time_cp, double **U_cp,
                        double **V_cp, double **TEMP_cp, int imax, int jmax)
{
    
    for(int i=0; i<imax;i++)
    {
        for(int j=0;j<jmax;j++)
        {
            U_cp[i][j] = U[i][j];
            V_cp[i][j] = V[i][j];
            TEMP_cp[i][j] = TEMP[i][j];
        }
    }
    
}


void restore_checkpoint(double *time, double **U, double **V, double **TEMP, double time_cp, double **U_cp,
                        double **V_cp,double **TEMP_cp, int imax, int jmax)
{
    
        for (int i=0;i<imax;i++)
        {
            for(int j=0; j<jmax;j++)
            {
                U[i][j] = U_cp[i][j];
                V[i][j] = V_cp[i][j];
                TEMP[i][j] = TEMP_cp[i][j];
            }
        }
    

}
*/
