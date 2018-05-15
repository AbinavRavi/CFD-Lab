#include <math.h>
#include "uvp.h"
#include"helper.h"
#include"boundary_val.h"
#include <stdio.h>

/////////////////////////////////////////////////////////////////////
void calculate_dt(double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,double Pr, int include_temp)
{

 double Umax = fabs(U[0][0]);
 
   for( int i = 0 ; i < imax ; i++ )
   {
      for( int j = 0 ; j < jmax ; j++ )
      {
         if ( fabs(U[i][j]) > Umax )
            Umax = fabs(U[i][j]);
      }
   }
   
  double Vmax = fabs(V[0][0]);
 
   for( int i = 0 ; i < imax ; i++ )
   {
      for( int j = 0 ; j < jmax ; j++ )
      {
         if ( fabs(V[i][j]) > Vmax )
            Vmax = fabs(V[i][j]);
      }
   }
 
 double x, y, z, w;
 
 x = (Re/2.0)*pow(((1.0/pow(dx,2.0))+(1.0/pow(dy,2.0))),-1.0);

 y = dx/(fabs(Umax));

 z = dy/(fabs(Vmax));
 
  double min = fmin(x,y);
    min = fmin(min,z);

 if(include_temp){
 w = ((Re*Pr)/2.0)*( (dx*dx)*(dy*dy)/( (dx*dx)+(dy*dy) ) );
    min = fmin(min,w);
 }

 
 if( (tau > 0) && (tau < 1))
 {
   *dt = tau * min;
 }
}

/////////////////////////////////////////////////////////////////////////////////
void calculate_fg(double Re,
		 double GX, double GY,
		 double alpha,
		 double dt,
		 double dx, double dy,
		 int imax, int jmax,
		 double** U, double** V,
		 double** F, double** G,int **flag,
		 double beta, double** temp, int include_temp)
{



 /*set boundary values in case of no-slip/free-slip/and inflow*/
for(int i = 0; i<imax; ++i)
{
	for(int j = 0; j<jmax; ++j)
	{
		if ( B_O(flag[i][j]) )  F[i][j] = U[i][j];

		if ( B_W(flag[i][j]) )  F[i-1][j] = U[i-1][j];

		if ( B_N(flag[i][j]) )  G[i][j] = V[i][j];

		if ( B_S(flag[i][j]) )  G[i][j-1] = V[i][j-1];

		if ( B_NO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j] = V[i][j]; }

		if ( B_NW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j] = V[i][j]; }

		if ( B_SO(flag[i][j]) ) { F[i][j] = U[i][j]; G[i][j-1] = V[i][j-1]; }

		if ( B_SW(flag[i][j]) ) { F[i-1][j] = U[i-1][j]; G[i][j-1] = V[i][j-1]; }
		//Outflow
		//if (flag[i][j]&(1<<3) ) P[i][j] = 0;
		//Inflow
		if (flag[i][j]&(1<<4) ) F[i][j] = U[i][j];
	}
}


    for(int i=0; i<imax-1; i++)
    {
        for(int j=0; j<jmax; j++)
	{
	if( ((flag[i][j]&(1<<0))&flag[i+1][j]) || ( (flag[i+1][j] & (1<<3)) && (flag[i][j]&(1<<0))) )
	//if((flag[i][j]&(1<<0))&flag[i+1][j]||((flag[i+1][j]&(1<<3)) && (flag[i][j]&(1<<0))))
	{
	if(include_temp)
	{
        F[i][j]=U[i][j]+dt*(
                //Central difference scheme for second derivatives
                (1/Re)*((U[i-1][j]-2*U[i][j]+U[i+1][j])/pow(dx,2.0)+(U[i][j-1]-2*U[i][j]+U[i][j+1])/pow(dy,2.0))
                //Modified Donor Cell method on convective term (d(u^2)/dx)
                -(1/dx)*0.25*(
                        (pow((U[i+1][j]+U[i][j]),2.0) - pow((U[i-1][j]+U[i][j]),2.0))
                        +alpha*(fabs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                             )
                 //Modified Donor Cell method on convective term (d(uv)/dy)
                -(1/dy)*0.25*(
                        ((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])- (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))
                    +alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))
                             )+GX)
                //Gravity component in x-direction
                -GX*beta*dt*(temp[i][j]+temp[i+1][j])/2;
	}
	else
	{
        F[i][j]=U[i][j]+dt*(
                //Central difference scheme for second derivatives
                (1/Re)*((U[i-1][j]-2*U[i][j]+U[i+1][j])/pow(dx,2.0)+(U[i][j-1]-2*U[i][j]+U[i][j+1])/pow(dy,2.0))
                //Modified Donor Cell method on convective term (d(u^2)/dx)
                -(1/dx)*0.25*(
                        (pow((U[i+1][j]+U[i][j]),2.0) - pow((U[i-1][j]+U[i][j]),2.0))
                        +alpha*(fabs(U[i+1][j]+U[i][j])*(U[i][j]-U[i+1][j])-fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                             )
                 //Modified Donor Cell method on convective term (d(uv)/dy)
                -(1/dy)*0.25*(
                        ((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])- (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))
                    +alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))
                             )+GX);		
	}

	}
	}
    }


    for(int i=0; i<imax; i++)
	{
        for(int j=0; j<jmax-1; j++)
	{
	if((flag[i][j]&(1<<0))&flag[i][j+1])
	{
	if(include_temp)
	{
        G[i][j]=V[i][j]+dt*(
                //Central difference Scheme for second derivatives
                (1/Re)*((V[i-1][j]-2*V[i][j]+ V[i+1][j])/pow(dx,2.0)+(V[i][j-1]-2*V[i][j]+ V[i][j+1])/pow(dy,2.0))
                //Modified Donor cell method for d(uv)/dx
                -(1/dx)*0.25*(
                        ((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])- (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))
                    +alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))
                            )
                //modified Donor cell method for d(v^2)/dy
                -(1/dy)*0.25*(
                        (pow((V[i][j]+V[i][j+1]),2.0) - pow((V[i][j-1]+V[i][j]),2.0))
                        +alpha*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
                            )+GY)
                -GY*beta*dt*(temp[i][j]+temp[i][j+1])*0.5;
	}
	else
	{
        G[i][j]=V[i][j]+dt*(
                //Central difference Scheme for second derivatives
                (1/Re)*((V[i-1][j]-2*V[i][j]+ V[i+1][j])/pow(dx,2.0)+(V[i][j-1]-2*V[i][j]+ V[i][j+1])/pow(dy,2.0))
                //Modified Donor cell method for d(uv)/dx
                -(1/dx)*0.25*(
                        ((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])- (U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))
                    +alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))
                            )
                //modified Donor cell method for d(v^2)/dy
                -(1/dy)*0.25*(
                        (pow((V[i][j]+V[i][j+1]),2.0) - pow((V[i][j-1]+V[i][j]),2.0))
                        +alpha*(fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
                            )+GY);
	
	}
	}
        }

    }
    //printf("F,G calculated \n");

}


////////////////////////////////////////////////////////////////////////////////
void calculate_uv(double dt,double dx,double dy,int imax, int jmax,
		 double**U, double**V,double**F,double**G,double **P,int **flag)
{
	for (int i = 0; i< imax-1;i++)
	{
		for (int j = 0; j<jmax;j++)
		{
			if(((flag[i][j]&(1<<0))&flag[i+1][j]) || ( (flag[i+1][j] & (1<<3)) && (flag[i][j]&(1<<0))))
			//update the U component of velocity
			U[i][j] = F[i][j] - (dt/dx)*(P[i+1][j]-P[i][j]);
		}
	}
	
	for (int i = 0; i< imax;i++)
	{
		for (int j = 0; j<jmax-1;j++)
		{
			if((flag[i][j]&(1<<0))&flag[i][j+1])
			//update the V component of velocity
			V[i][j] = G[i][j] - (dt/dy)*(P[i][j+1] -P[i][j]);
		}
	}
    //printf("U,V calculated \n");

}

/////////////////////////////////////////////////////////////////
void calculate_rs(double dt,
		  double dx,
		  double dy,
		  int imax,
		  int jmax,
		  double **F,
		  double **G,
		  double **RS,int **flag)
{

	for(int i=0; i<imax; i++)
   	{
        	for(int j=0; j<jmax; j++)
		{
			if(flag[i][j]&(1<<0))
			RS[i][j] =  (1/dt)*( (F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy );
		}
	}
    //printf("RS calculated \n");
}

void calculate_temp(double **temp, double **temp1, double Pr, double Re, int imax,int jmax,double dx, double dy,
		double dt, double alpha,double **U,double **V,int **flag, double TI, double T_h, double T_c, int select)
{   
    for(int i = 0; i<imax; ++i)
    {
  	for(int j = 0; j<jmax; ++j)
  	{
  		if ( B_O(flag[i][j]) )  temp[i][j] = temp[i+1][j];

  		if ( B_W(flag[i][j]) )  temp[i][j] = temp[i-1][j];

  		if ( B_N(flag[i][j]) )  temp[i][j] = temp[i][j+1];

  		if ( B_S(flag[i][j]) )  temp[i][j] = temp[i][j-1];

  		if ( B_NO(flag[i][j]) ) temp[i][j] = (temp[i][j+1] + temp[i+1][j])/2;

  		if ( B_NW(flag[i][j]) ) temp[i][j] = (temp[i][j+1] + temp[i-1][j])/2;

  		if ( B_SO(flag[i][j]) ) temp[i][j] = (temp[i][j-1] + temp[i+1][j])/2;

  		if ( B_SW(flag[i][j]) ) temp[i][j] = (temp[i][j-1] + temp[i-1][j])/2;

  		if (flag[i][j]&(1<<3) ) temp[i][j] = temp[i-1][j];

  		if (flag[i][j]&(1<<4) ) temp[i][j] = TI;
  	}
  }

switch(select)
{
	case 3:
	for(int j=0; j<jmax; j++)
	{
		temp[0][j] = 2*T_h - temp[1][j];
		temp[imax-1][j] = 2*T_c - temp[imax-2][j];		
	}
	break;
	
	case 4:
	for(int j=0; j<jmax; j++)
	{
		temp[0][j] = 2*T_h - temp[1][j];
		temp[imax-1][j] = 2*T_c - temp[imax-2][j];		
	}
	break;

	case 5:
	for(int i=0; i<imax; i++)
	{
		temp[i][0] = 2*T_h - temp[i][1];
		temp[i][jmax-1] = 2*T_c - temp[i][jmax-2];		
	}

}

double dut_dx;
double dvt_dy;
double dt2_dx2;
double dt2_dy2;
double Z;
for (int i = 0; i< imax;i++)
{
    for (int j=0; j<jmax;j++)
	{

	if(flag[i][j]&((1<<0)|(1<<3)|(1<<4)))
	{

  		dut_dx = (1/dx)*( (U[i][j]*(temp[i][j]+temp[i+1][j])*0.5)-(U[i-1][j]*(temp[i-1][j]+temp[i][j])*0.5))+(alpha/dx)*(
				(fabs(U[i][j])*(temp[i][j]-temp[i+1][j])*0.5) - (fabs(U[i-1][j])*(temp[i-1][j]-temp[i][j])*0.5)
				);

  		dvt_dy = (1/dy)*( (V[i][j]*(temp[i][j]+temp[i][j+1])*0.5)-(V[i][j-1]*(temp[i][j-1]+temp[i][j])*0.5) )+(alpha/dy)*(
				(fabs(V[i][j])*(temp[i][j]-temp[i][j+1])*0.5) - (fabs(V[i][j-1])*(temp[i][j-1]-temp[i][j])*0.5)
				);

  		dt2_dx2 = (temp[i+1][j] - 2*temp[i][j] + temp[i-1][j])/(dx*dx);

  		dt2_dy2 = (temp[i][j+1] - 2*temp[i][j] +temp[i][j-1])/(dy*dy);

  		Z = (1/(Re*Pr))*(dt2_dx2+dt2_dy2) - dut_dx - dvt_dy;

      	temp1[i][j] = temp[i][j]+ (dt*Z);
    }
  }
	
}

  for (int i = 0; i< imax;i++){
    for (int j=0;j<jmax;j++){

		if(flag[i][j]&((1<<0)|(1<<3)|(1<<4))){

		temp[i][j] = temp1[i][j];
			//printf("%d ",(int)temp[i][j]);
    }
  }
//printf("\n");	
}

}

///////////////////////////////////////////////////////////////////////////////

void nullify_obstacles1(double **U, double **V, double **P, int **flag, int imax, int jmax)
{

	for (int i = 0; i< imax;i++)
	{
	    for (int j=0;j<jmax; j++)
		{
			if(flag[i][j]&( (1<<1)|(1<<2)) )
			{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;
			}
			
		}
	}

}


void nullify_obstacles2(double **U, double **V, double **P, double **T, int **flag, int imax, int jmax)
{
	for (int i = 0; i< imax;i++)
	{
	    for (int j=0;j<jmax; j++)
		{
			if(flag[i][j]&( (1<<1)|(1<<2)) )
			{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;
				T[i][j] = 0;
			}
		}
	}


}

