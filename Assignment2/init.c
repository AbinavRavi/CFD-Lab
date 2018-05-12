#include "helper.h"
#include "init.h"
#include <stdio.h>

int read_parameters( const char *szFileName,       /* name of the file */
		    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/ 
		    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
		    double *dt,                /* time step */
		    double *t_end,             /* end time */
		    double *tau,               /* safety factor for time step*/
		    double *dt_value,		/* time for output */
		    double *eps,               /* accuracy bound for pressure*/
		    double *omg,               /* relaxation factor */
		    double *alpha,             /* uppwind differencing factor*/
                    int  *itermax,             /* max. number of iterations  */
		    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
		    double *Re,                /* reynolds number   */
                    double *Pr,
		    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
       		    double *TI,
		    double *T_h,
		    double *T_c,
		    double *beta,
		    double *dx,                /* length of a cell x-dir. */
                    double *dy                /* length of a cell y-dir. */
                    		                               /* for pressure per time step */

)           
{
   printf("PROGRESS: Reading .dat file... \n");
   //READ_STRING( szFileName, *problem );
   //READ_STRING( szFileName, geometry );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *dt    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *alpha );
   READ_INT   ( szFileName, *itermax );

   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *Re );
   READ_DOUBLE( szFileName, *Pr );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *TI );
   READ_DOUBLE( szFileName, *T_h );
   READ_DOUBLE( szFileName, *T_c );
   READ_DOUBLE( szFileName, *beta );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   printf("PROGRESS: .dat file read... \n \n");
 return 1;

}


void init_uvp(double UI, double VI, double PI, int imax, int jmax,
		 double** U, double** V, double** P, int** flag)
{
	printf("PROGRESS: Starting initialization of U,V,P ... \n");
	
	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if(flag[i][j]&(1<<0)){

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
			}
			/*else{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;		
			}*/
		}
	}
	printf("PROGRESS: U,V,P matrices initialized... \n \n");
}

void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax,
		 double** U, double** V, double** P, double** T, int** flag)
{
	printf("PROGRESS: Starting initialization of U,V,P,T ... \n");
	
	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if(flag[i][j]&(1<<0)){

				U[i][j] = UI;
				V[i][j] = VI;
				P[i][j] = PI;
				T[i][j] = TI;
			}
			/*else{
				U[i][j] = 0;
				V[i][j] = 0;
				P[i][j] = 0;
				T[i][j] = 0;		
			}*/
		}
	}
	printf("PROGRESS: U,V,P,T matrices initialized... \n \n");
}

int  isfluid(int pic){
	if((pic == 2)||(pic == 3)||(pic == 4)) {return 1;}
		else {return 0;}
}

//int assert(int imax, int jmax, int **pic){


//}


void init_flag(char* problem, char* geometry, int imax, int jmax, int **flag)
{
	printf("PROGRESS: Setting flags... \n");
	int **pic = read_pgm(geometry);

	for (int i=0; i<imax; i++)
	{
		for (int j=0; j<jmax; j++)
		{

		flag[i][j] = 0;

		//assert(pic);
		switch(pic[i][j])
		{
			case 0: //fluid
			flag[i][j] = 1<<1;
			break;

			case 1: //no-slip
			flag[i][j] = 1<<2;
			break;

			case 2: //free-slip
			flag[i][j] = 1<<3;
			break;

			case 3: //outflow
			flag[i][j] = 1<<4;
			break;

			case 4: //inflow
			flag[i][j] = 1<<0;
			break;
		}

			if(!isfluid(pic[i][j])) //set boundaries if not 
			{

				if(i<imax-1 && pic[i+1][j]==4)
				{
				flag[i][j] |= 1<<8;
				}
				if( i>0 && pic[i-1][j]==4)
				{
				flag[i][j] |= 1<<7;
				}
				if(j<jmax-1 && pic[i][j+1]==4)
				{
				flag[i][j] |= 1<<5;
				}
				if(j>0 && pic[i][j-1]==4)
				{
				flag[i][j] |= 1<<6;
				}
			}
 

		}

	}
	/*printf("\n");
	for(int j=0;j<jmax;j++){

		for(int i=0;i<imax;i++){

			printf("%d",pic[i][j]);
		}
		printf("\n");	
	}
	for(int j=0;j<jmax;j++){

		for(int i=0;i<imax;i++){

			printf("%d ",flag[i][j]);
		}
		printf("\n");	
	}*/
	printf("PROGRESS: flags set using .pgm file...\n \n");


}
