

struct particle{
    double x, y;
    struct particle *next;
};

struct particleline{
    int length;
    struct particle *Particles;
};

struct particleline *Particlelines;


struct particleline *INIT_PARTICLES (int *N, int imax, int jmax, double delx, double dely, int ppc, int **flag){
    
    *N = 0;
    int len = sqrt(ppc);
    for(int i = 0; i < imax; i++)
    {
        for(int j = 0; j < jmax; j++)
        {
            if (flag[i][j]&&(1<<0|1<<3|1<<4)!=0) {
                
                for(int k = 0; k < len; k++)
                {
                    (*N)++;
                    Particlelines[(*N)+k].length = len;
                    
                    for(int p = 0; p < len; p++)
                    {
                        Particlelines[(*N)+k].Particles[p].x = (i-1)*delx + (delx/len)*p;
                        Particlelines[(*N)+k].Particles[p].y = (j-1)*dely + (dely/len)*k;
                        Particlelines[(*N)+k].Particles[p].next = (p<len-1)?&Particlelines[(*N)+k].Particles[p+1]:0;
                    }

                }   
                
            }
        }
        
    }

    return &Particlelines[0];
}
    


void MARK_CELLS(int **flag, int imax, int jmax, double delx, double dely, int N, struct particleline *Partlines){

    for(int i = 0; i < imax; i++)
    {
        for(int j = 0; j < jmax; j++)
        {
            if (flag[i][j]&&(1<<1|1<<2)==0) {
                
                for(int k = 0; k < N; k++)
                {   
                    for(int p = 0; p < Particlelines[k].length; p++)
                    {
                        //Set flag for inflow, outflow, fluid or empty region
                        if( (Particlelines[k].Particles[p].x >= (i-1)*delx)&&(Particlelines[k].Particles[p].x< i*delx) 
                        && (Particlelines[k].Particles[p].y >= (j-1)*dely)&&(Particlelines[k].Particles[p].y< j*dely) )
                        {
                            if (i == 0)             flag[i][j] = 1<<4;
                            else if (i == imax-1)   flag[i][j] = 1<<3;
                            else                    flag[i][j] = 1<<0;
                        }
                        else
                        {
                            flag[i][j] = 1<<11;
                        }
                        //set flag for surface or interior region, set free boundaries
                        if( (Particlelines[k].Particles[p].x >= (i-1)*delx)&&(Particlelines[k].Particles[p].x< i*delx) 
                        && (Particlelines[k].Particles[p].y >= (j-1)*dely)&&(Particlelines[k].Particles[p].y< j*dely) )
                        {
                            if (i<imax-1 && (flag[i+1][j]&&(1<<11)!=0) )
                            {
                                flag[i][j] = 1<<14;
                                flag[i][j] |= 1<<10;
                            }
                            if (i>0 && (flag[i-1][j]&&(1<<11)!=0) )
                            {
                                flag[i][j] = 1<<15;
                                flag[i][j] |= 1<<10;
                            }
                            if (j<jmax-1 && (flag[i][j+1]&&(1<<11)!=0) )
                            {
                                flag[i][j] = 1<<13;
                                flag[i][j] |= 1<<10;
                            }
                            if (j>0 && (flag[i][j-1]&&(1<<11)!=0) )
                            {
                                flag[i][j] = 1<<12;
                                flag[i][j] |= 1<<10;
                            }

                            if (flag[i][j]&&(1<<10)==0) {
                                flag[i][j] |= 1<<9;
                            }
                        
                         }

                         //Set obstacle boundaries
                         if (flag[i][j]&&(1<<1|1<<2)!=0)
                         {
                            if (i<imax-1 && (flag[i+1][j]&&(1<<0|1<<3|1<<4)!=0) )
                            {
                                flag[i][j] = 1<<8;
                            }
                            if (i>0 && (flag[i-1][j]&&(1<<0|1<<3|1<<4)!=0) )
                            {
                                flag[i][j] = 1<<7;
                            }
                            if (j<jmax-1 && (flag[i][j+1]&&(1<<0|1<<3|1<<4)!=0) )
                            {
                                flag[i][j] = 1<<5;
                            }
                            if (j>0 && (flag[i][j-1]&&(1<<0|1<<3|1<<4)!=0) )
                            {
                                flag[i][j] = 1<<6;
                            }
                         }
                         
                    }

                }   
                
            }
        }
        
    }


}

void SET_UVP_SURFACE(double **U,double **V, double **P, int **flag, int imax, int jmax, double Re, double delx, double dely, double delt)
{
    
}