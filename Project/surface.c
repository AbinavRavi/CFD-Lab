#include "surface.h"


struct particleline *INIT_PARTICLES (int *N, int imax, int jmax, double delx, double dely, int ppc, int **flag){
    
    *N = 0;
    int len = sqrt(ppc);

    // find the total number of particlelines
    for(int i = 0; i < imax; i++)
    {
        for(int j = 0; j < jmax; j++)
        {
            if ((flag[i][j]&(1<<0|1<<3|1<<4))!=0) 
            {    
                    (*N)+= len;

            }
        }
    }
    printf("Number of particlelines: %d \n",*N);

    //allocate particleline to an array of size N
    struct particleline *Particlelines = (struct particleline*)malloc(sizeof(struct particleline)*(*N));
    
    int n = 0;
    for(int i = 0; i < imax; i++)
    {
        for(int j = 0; j < jmax; j++)
        {
            
            if ((flag[i][j]&(1<<0|1<<3|1<<4))!=0)
            {    
                for(int k = 0; k < len; k++)
                {
                    struct particle *start = (struct particle *)malloc(sizeof(struct particle));
                    
                    Particlelines[n + k].length = len;
                    Particlelines[n + k].Particles = start;
                    insert_particles(start, i, j, delx, dely, len, k);
                }
                n += len;
                
            }
            
        } 
    }

    return &Particlelines[0];
}

void insert_particles(struct particle *start, int i, int j, double delx, double dely, int len, int cell_count){

    struct particle *list;

    list = start;
    start->next = NULL;

    for(int k=0; k<len; k++)
    {   
        list->x = i*delx + (delx/(double)len)*k;
        list->y = j*dely + (dely/(double)len)*cell_count;
        
        if (k<len-1) {
            list->next = (struct particle *)malloc(sizeof(struct particle));
        }
        else  list->next = NULL;
        
        list = list->next;
    }   
    
   
}

void MARK_CELLS(int **flag, int imax, int jmax, double delx, double dely, int N, struct particleline *Partlines){

    for(int i = 0; i < imax; i++)
    {
        for(int j = 0; j < jmax; j++)
        {
            if ( (flag[i][j]&((1<<1|1<<2)))==0 )//not obstacles
            {
                flag[i][j] = 1<<11; // preset all open cells to empty

                for(int k = 0; k < N; k++)
                {
                    struct particle *temp = Partlines[k].Particles;

                    for(int p = 0; p < Partlines[k].length; p++)
                    {
                        //Set flag for inflow, outflow, fluid or empty region
                        if( ((temp->x) >= i*delx)&&((temp->x) < (i+1)*delx) 
                        &&  ((temp->y) >= j*dely)&&((temp->y) < (j+1)*dely)  )
                        {
                            if (i == 0)             {flag[i][j] = 1<<4;} //inflow
                            else if (i == imax-1)   {flag[i][j] = 1<<3;} //outflow
                            else                    {flag[i][j] = 1<<0;} //fluid
                        }
                        if( ((temp->x) >= i*delx)&&((temp->x) < (i+1)*delx) 
                        &&  ((temp->y) >= j*dely)&&((temp->y) < (j+1)*dely) )
                        {
                            //set flag for surface or interior region, set free boundaries
                            if (i<imax-1 && ((flag[i+1][j]&(1<<11))!=0) )
                            {
                                flag[i][j] |= 1<<15; //east
                                flag[i][j] |= 1<<10;
                            }
                            if (i>0 && ((flag[i-1][j]&(1<<11))!=0) )
                            {
                                flag[i][j] |= 1<<14; //west
                                flag[i][j] |= 1<<10;
                            }
                            if (j<jmax-1 && ((flag[i][j+1]&(1<<11))!=0) )
                            {
                                flag[i][j] |= 1<<12; //north
                                flag[i][j] |= 1<<10;
                            }
                            if (j>0 && ((flag[i][j-1]&(1<<11))!=0) )
                            {
                                flag[i][j] |= 1<<13; //south
                                flag[i][j] |= 1<<10;
                            }

                            if ((flag[i][j]&(1<<10))==0) { //interior cells
                                flag[i][j] |= 1<<9;
                            }
                        
                         }
                        //printf("Debug1 \n");
                         temp = temp->next;
                         //printf("Debug2 \n");
                    }
                }
            }
            //Set obstacle boundaries
            else
            {  
                if (i<imax-1 && (flag[i+1][j]&(1<<0|1<<3|1<<4)) )
                {
                    flag[i][j] |= 1<<8; //east
                }
                if (i>0 && (flag[i-1][j]&(1<<0|1<<3|1<<4)) )
                {
                    flag[i][j] |= 1<<7; //west
                }
                if (j<jmax-1 && (flag[i][j+1]&(1<<0|1<<3|1<<4)) )
                {
                    flag[i][j] |= 1<<5; //north
                }
                if (j>0 && (flag[i][j-1]&(1<<0|1<<3|1<<4)) )
                {
                    flag[i][j] |= 1<<6; //south
                }
            }


        }  
    }


}

//Surface type 1
int S_O(int flag){
	return ( (flag & (1<<15)) && ~( flag & ((1<<12)|(1<<13)|(1<<14)) ) );
}

int S_W(int flag){
	return ( (flag & (1<<14)) && ~( flag & ((1<<12)|(1<<13)|(1<<15)) ) );
}

int S_N(int flag){
	return ( (flag & (1<<12)) && ~( flag & ((1<<15)|(1<<13)|(1<<14)) ) );
}

int S_S(int flag){
	return ( (flag & (1<<13)) && ~( flag & ((1<<12)|(1<<15)|(1<<14)) ) );
}

// Surface type 2
int S_NO(int flag){
	return ( (flag & (1<<12)) && (flag & (1<<15)) &&  ~(flag & ((1<<13)|(1<<14))) );
}

int S_OS(int flag){
	return ( (flag & (1<<15)) && (flag & (1<<13)) &&  ~(flag & ((1<<12)|(1<<14))) );
}

int S_SW(int flag){
	return ( (flag & (1<<13)) && (flag & (1<<14)) &&  ~(flag & ((1<<12)|(1<<15))) );
}

int S_WN(int flag){
	return ( (flag & (1<<14)) && (flag & (1<<12)) &&  ~(flag & ((1<<13)|(1<<15))) );
}

//Surface type 3
int S_OW(int flag){
	return ( (flag & (1<<14)) && (flag & (1<<15)) &&  ~(flag & ((1<<13)|(1<<12))) );
}

int S_NS(int flag){
	return ( (flag & (1<<12)) && (flag & (1<<13)) &&  ~(flag & ((1<<14)|(1<<15))) );
}

//Surface type 4
int S_NOS(int flag){
	return ( (flag & (1<<12)) && (flag & (1<<15)) && (flag & (1<<13)) && ~(flag & (1<<14)) );
}

int S_OSW(int flag){
	return ( (flag & (1<<15)) && (flag & (1<<13)) && (flag & (1<<14)) && ~(flag & (1<<12)) );
}

int S_SWN(int flag){
	return ( (flag & (1<<13)) && (flag & (1<<14)) && (flag & (1<<12)) && ~(flag & (1<<15)) );
}

int S_WNO(int flag){
	return ( (flag & (1<<14)) && (flag & (1<<12)) && (flag & (1<<15)) && ~(flag & (1<<13)) );
}

//Surface type 5
int S_NOSW(int flag){
	return ( (flag & (1<<14)) && (flag & (1<<12)) && (flag & (1<<15)) && (flag & (1<<13)) );
}



void SET_UVP_SURFACE(double **U,double **V, double **P, int **flag, int imax, int jmax, double Re, double delx, double dely, double delt,double GX, double GY)
{
    for (int i=0; i< imax; i++)
    {
        for (int j=0;j<jmax;j++)
        {
    //surface type 1
            if (S_O(flag[i][j]))
            {
                //printf("DebugO \n");
                U[i][j] = U[i-1][j] - ((delx/dely)*(V[i][j]-V[i][j-1]));
                P[i][j] = (2/Re)*(U[i][j]-U[i-1][j])/delx;
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                    V[i+1][j-1] = V[i][j-1] -((delx/dely)*(U[i][j] - U[i][j-1]));
                }

            }
            if (S_W(flag[i][j]))
            {
                //printf("DebugW \n");
                U[i-1][j] = U[i][j] +(delx/dely)*(V[i][j] - V[i][j-1]);
                P[i][j] = (2/Re)*(U[i][j] - U[i-1][j]);
                if((flag[i-1][j-1]&(1<<11))!=0)
                {   
                    V[i-1][j-1] = V[i][j-1]+(delx/dely)*(U[i-1][j] - U[i-1][j-1]);
                }

            }
            if (S_N(flag[i][j]))
            {
               // printf("DebugN (%d %d) \n", i, j);
                V[i][j] = V[i][j-1] - (dely/delx)*(U[i][j]-U[i-1][j]);
                //printf("Debug \n");
                P[i][j] = (2/Re)*((V[i][j]-V[i][j-1])/dely);
               // printf("Debug \n");
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    //printf("Debug \n");
                    U[i-1][j+1] = U[i][j] - (dely/delx)*(V[i][j] - V[i][j-1]);
                }
            }
            if (S_S(flag[i][j]))
            {
               // printf("DebugS \n");
                V[i][j-1] = V[i][j]+ (dely/delx)*(U[i][j] - U[i-1][j]);
                P[i][j] = (2/Re)*((V[i][j] - V[i][j-1])/dely);
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j] + (dely/delx)*(V[i][j-1]-V[i-1][j-1]);
                }
            }
            //surface type 2

            if (S_NO(flag[i][j]))
            {
               // printf("DebugNO \n");
                U[i][j] = U[i-1][j];
                V[i][j] = V[i][j-1];
                P[i][j] = (0.5/Re)*(((U[i-1][j]+U[i][j]-U[i-1][j-1] - U[i][j-1])/dely)+((V[i][j-1]+V[i][j]-V[i-1][j]-V[i-1][j-1])/delx));
                if ((flag[i+1][j-1]&(1<<11))!=0)
                {
                    U[i][j-1] = U[i][j] +(delx/dely)*(V[i+1][j-1]-V[i][j-1]);
                }
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    V[i-1][j] = V[i][j] +(dely/delx)*(U[i-1][j+1]-U[i-1][j]);
                }
                U[i][j+1] = U[i][j];
                V[i+1][j] = V[i][j];

            }

            if(S_OS(flag[i][j]))
            {
                U[i][j] = U[i-1][j];
                V[i][j-1] = V[i][j];
                P[i][j] = (-0.5/Re)*(((U[i][j+1]+U[i-1][j+1]-U[i][j]-U[i-1][j])/dely)+((V[i][j]+V[i][j-1]-V[i-1][j]-V[i-1][j-1])/delx));
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j] +(dely/delx)*(V[i][j-1] - V[i-1][j-1]);
                }
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                    V[i+1][j] = V[i][j] - (delx/dely)*(U[i][j+1] - U[i][j]);
                }
                U[i][j-1] = U[i][j];
                V[i+1][j-1] = V[i][j-1];
            }   
            if(S_SW(flag[i][j]))
            {
                U[i-1][j] = U[i][j];
                V[i][j-1] = V[i][j];
                P[i][j] = (0.5/Re)*(((U[i][j+1]+U[i-1][j+1]-U[i][j] - U[i-1][j])/dely)+((V[i][j]+V[i][j-1]-V[i+1][j]-V[i+1][j-1])/delx));
                if ((flag[i+1][j-1]&(1<<11))!=0)
                {
                    U[i][j-1] = U[i][j] +(dely/delx)*(V[i+1][j-1]-V[i][j-1]);
                }
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    V[i-1][j] = V[i][j] +(delx/dely)*(U[i-1][j+1]-U[i-1][j]);
                }
                U[i-1][j-1] = U[i-1][j];
                V[i-1][j-1] = V[i][j-1];

            }
            if(S_WN(flag[i][j]))
            {
                U[i-1][j] = U[i][j];
                V[i][j] = V[i][j-1];
                P[i][j] = (-0.5/Re)*(((U[i-1][j]+U[i][j]-U[i][j-1]-U[i-1][j-1])/dely)+((V[i+1][j]+V[i+1][j-1]-V[i][j]-V[i][j-1])/delx));
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j]+(delx/dely)*(V[i][j-1]-V[i-1][j-1]);
                }
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                    V[i+1][j] = V[i][j] - (dely/delx)*(U[i][j+1] - U[i][j]);
                }
                U[i-1][j+1] = U[i-1][j];
                V[i+1][j] = V[i][j];

            }

            //surface type 3
            if(S_OW(flag[i][j]))
            {
                U[i][j] = U[i][j]+(delt*GX);
                U[i-1][j] = U[i-1][j] +(delt*GX);
                P[i][j] = 0;
                if(((flag[i][j]&(1<<11))!=0))
                {
                    V[i-1][j-1] = V[i][j-1] +(delx/dely)*(U[i-1][j] - U[i-1][j-1]);
                }
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                    V[i+1][j-1] = V[i][j-1] - (delx/dely)*(U[i][j] - U[i][j-1]);
                }

            }

            if (S_NS(flag[i][j]))
            {
                V[i][j] = V[i][j]+(delt*GY);
                V[i][j-1] = V[i][j]+(delt*GY);
                P[i][j] = 0;
                if ((flag[i-1][j+1]&(1<<11))!=0)
                {
                    U[i-1][j+1] = U[i-1][j] - (dely/delx)*(V[i][j] - V[i-1][j]);
                }
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j]+(dely/delx)*(V[i][j-1]-V[i-1][j-1]);
                }

            }

            //surface type 4
            if(S_NOS(flag[i][j]))
            {
                P[i][j] = 0;
                V[i][j] = V[i][j]+(delt*GY);
                V[i][j-1] = V[i][j-1]+(delt*GY);
                U[i][j] = U[i-1][j] -(delx/dely)*(V[i][j] - V[i][j-1]);
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                    U[i][j+1] = U[i][j];
                    V[i+1][j] = V[i][j];
                }
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                    U[i][j-1] = U[i][j];
                    V[i+1][j-1] = V[i][j-1];
                }
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j+1] = U[i-1][j]- (dely/delx)*(V[i][j] - V[i-1][j]);
                }
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                     U[i-1][j-1] = U[i-1][j]+ (dely/delx)*(V[i][j-1] - V[i-1][j-1]);
                }

            }
            if(S_OSW(flag[i][j]))
            {
                P[i][j] = 0;
                U[i][j] = U[i][j]+(delt*GX);
                U[i-1][j] = U[i-1][j] +(delt*GX);
                V[i][j-1] = V[i][j] - (dely/delx)*(U[i][j]-U[i-1][j]);
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j];
                    V[i-1][j-1] = V[i][j-1];
                }
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                    U[i][j-1] = U[i][j];
                    V[i+1][j-1] = V[i][j-1];
                }
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    V[i-1][j] = V[i][j]+(delx/dely)*(U[i-1][j+1]-U[i-1][j]);
                }
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                    V[i+1][j] = V[i][j] - (delx/dely)*(U[i][j+1] - U[i][j]);
                }

            }
            if(S_SWN(flag[i][j]))
            {
                P[i][j] = 0;
                V[i][j] = V[i][j]+(delt*GY);
                V[i][j-1] = V[i][j-1]+(delt*GY);
                U[i-1][j] = U[i][j] +(delx/dely)*(V[i][j] - V[i][j-1]);
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j];
                    V[i-1][j-1] = V[i][j-1];
                }
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    U[i-1][j+1] = U[i-1][j];
                    V[i-1][j] = V[i][j];
                }
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                    U[i][j-1] = U[i][j]+(dely/delx)*(V[i+1][j-1] - V[i][j-1]);
                }
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                     U[i][j+1] = U[i][j]- (dely/delx)*(V[i+1][j] - V[i][j]);
                }
            }
            if(S_WNO(flag[i][j]))
            {
                P[i][j] = 0;
                U[i-1][j] = U[i-1][j]+(delt*GX);
                U[i][j] = U[i][j] +(delt*GX);
                V[i][j] = V[i][j-1] - (dely/delx)*(U[i][j]-U[i-1][j]);
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    U[i-1][j+1] = U[i-1][j];
                    V[i-1][j] = V[i][j];
                }
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                    U[i][j+1] = U[i][j];
                    V[i+1][j] = V[i][j];
                }
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    V[i-1][j-1] = V[i][j-1]+(delx/dely)*(U[i-1][j]-U[i-1][j-1]);
                }
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                    V[i+1][j-1] = V[i][j-1] - (delx/dely)*(U[i][j] - U[i][j-1]);
                }

            }

            //surface type 5
            if(S_NOSW(flag[i][j]))
            {
                P[i][j] = 0;
                V[i][j] = V[i][j]+(delt*GY);
                V[i][j-1] = V[i][j-1]+(delt*GY);
                U[i][j] = U[i][j] +(delt*GX);
                U[i-1][j] = U[i-1][j]+(delt*GX);
                if((flag[i-1][j-1]&(1<<11))!=0)
                {
                    U[i-1][j-1] = U[i-1][j];
                    V[i-1][j-1] = V[i][j-1];
                }
                if((flag[i-1][j+1]&(1<<11))!=0)
                {
                    U[i-1][j+1] = U[i-1][j];
                    V[i-1][j] = V[i][j];
                }
                if((flag[i+1][j+1]&(1<<11))!=0)
                {
                   U[i][j+1] = U[i][j];
                   V[i+1][j] = V[i][j];
                }
                if((flag[i+1][j-1]&(1<<11))!=0)
                {
                   U[i][j-1] = U[i][j];
                   V[i+1][j-1] = V[i][j-1];
                }
                
            }
        }
    }

}

/*
void DELETE_PARTICLES(int **flag, int imax, int jmax, double delx, double dely, int N, struct particleline *Partlines){

    for(int i = 0; i < imax; i++)
    {
        for(int j = 0; j < jmax; j++)
        {
            if (flag[i][j]&(1<<0|1<<3|1<<4|1<<11)==0)
            {

                for(int k = 0; k < N; k++)
                {
                    struct particle *currP = NULL;
                    struct particle *prevP = NULL;

                    for(currP = Particlelines[k].Particles; currP != NULL; prevP = currP, currP = currP->next)
                    {
                        if( (currP->x >= i*delx)&&(currP->x < (i+1)*delx) 
                        &&  (currP->y >= j*dely)&&`(currP->y < (j+1)*dely) )
                        {
                           //unlink first particle 
                            if (currP == Particlelines[k].Particles) {
                                Particlelines[k].Particles = currP->next;
                            }
                            //unlink last particle
                            else if (currP->next == NULL) {
                                prevP->next = NULL;
                            }
                            //delete interior particle
                            else{
                                prevP->next = currP->next;
                            }
                            free(currP);
                            (Particlelines[k].length)--;
                        }
                    }
                }
            }
        }
    }
    
}
*/

double U_interp(double **U, double dx, double dy, double x, double y){

    double x1; double x2;
    double y1; double y2;

    double u1; double u2;
    double u3; double u4;
    
    int i = (int)(x/dx);
    int j = (int)((y+0.5*dy)/dy);

    x1 = i*dx;
    x2 = (i+1)*dx;
    y1 = (j-0.5)*dy;
    y2 = (j+0.5)*dy;

    
    if (i==0) { //inflow  

        u1 = 1;
        u2 = U[i][j-1];
        u3 = 1;
        u4 = U[i][j];
    }
    else{ //fluid or outflow

        u1 = U[i-1][j-1];
        u2 = U[i][j-1];
        u3 = U[i-1][j];
        u4 = U[i][j];
    }

    return (1/(dx*dy))*( (x2-x)*(y2-y)*u1 + (x-x1)*(y2-y)*u2
            +(x2-x)*(y-y1)*u3 + (x-x1)*(y-y1)*u4 );
}

double V_interp(double **V, double dx, double dy, double x, double y){
    
    double x1; double x2;
    double y1; double y2;

    double v1; double v2;
    double v3; double v4;

    int i = (int)((x+0.5*dx)/dx);
    int j = (int)(y/dy);
    
    x1 = (i-0.5)*dx;
    x2 = (i+0.5)*dx;
    y1 = j*dy;
    y2 = (j+0.5)*dy;


    if (i==0) { //inflow  

        v1 = 0;
        v2 = V[i][j-1];
        v3 = 0;
        v4 = V[i][j];
    }
    else{ //fluid or outflow

        v1 = V[i-1][j-1];
        v2 = V[i][j-1];
        v3 = V[i-1][j];
        v4 = V[i][j];
    }


    return (1/(dx*dy))*( (x2-x)*(y2-y)*v1 + (x-x1)*(y2-y)*v2
            +(x2-x)*(y-y1)*v3 + (x-x1)*(y-y1)*v4 );
}


void ADVANCE_PARTICLES(double **U, double **V, double delx, double dely, double delt, int N, struct particleline *Partlines, int **flag, int imax, int jmax){
    
    int i, j;
    for(int k = 0; k < N; k++)
    {
        
        if (Partlines[k].Particles != NULL)
        {
        struct particle *temp = Partlines[k].Particles;
        
        for(int p = 0; p < Partlines[k].length; p++)
        {
            double x_pos = temp->x;
            double y_pos = temp->y;
            
            i = (int)(x_pos/delx);
            j = (int)((y_pos+0.5*dely)/dely);

            if ((flag[i][j]&(1<<0|1<<3|1<<4))!=0) //only advance particles in fluid cells
            {
                
                if ((i>0)&&(i<imax)&&(j>0)&&(j<jmax)) {
                    temp->x += delt*U_interp(U, delx, dely, x_pos, y_pos);
                }
                
                //printf("%f \n", U_interp(U, delx, dely, x_pos, y_pos));
            }

            i = (int)((x_pos+0.5*delx)/delx);
            j = (int)(y_pos/dely);

            if ((flag[i][j]&(1<<0|1<<3|1<<4))!=0) //only advance particles in fluid cells
            {
                if ((i>0)&&(i<imax)&&(j>0)&&(j<jmax)) {
                    temp->y += delt*V_interp(V, delx, dely, x_pos, y_pos);
                }
            }
            //printf("(%f %f) \n",temp->x,temp->y);
            temp = temp->next;
        }
        }
    }

}


void FREE_PARTICLELINES(struct particleline *Partlines, int N, int imax, int jmax, double delx, double dely, int ppc, int **flag )
{
    for(int k = 0; k < N; k++)
    {
        struct particle *start = Partlines[k].Particles;
        free_particles(start);
    }
    free(Partlines);
}


void free_particles(struct particle *start)
{
   struct particle* current = start;
   struct particle* next;
 
   while (current != NULL) 
   {
       next = current->next;
       free(current);
       current = next;
   }
   
   start = NULL;
}