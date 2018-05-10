#include "helper.h"
#include "init.h"
#include <stdio.h>

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		    double *dt_value,		/* time for output */
		    char *problem,
                    char *geometry
                  double **temp
                double *Pr
                double *TI
                double *beta)           
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   READ_STRING( szFileName, problem );
   READ_STRING( szFileName, geometry );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   printf("Parameter file read \n");

   return 1;
}


void init_uvp(double UI, double VI, double PI, int imax, int jmax, double** U, double** V, double** P){
	init_matrix(U ,0,imax ,0, jmax, UI);
	init_matrix(V ,0,imax ,0, jmax, VI);
	init_matrix(P ,0,imax ,0, jmax, PI);
	printf("U,V,P matrices iniialized \n");
}

int  isfluid(int pic){
	if((pic == 2)||(pic == 3)||(pic == 4))
		return 1;
		else return 0;
}

void init_flag(const char* problem, const char* geometry, int imax, int jmax, int **flag)
{
	int **pic = read_pgm(geometry);

	for (int i=0; i<imax; i++)
	{
		for (int j=0; j<jmax; j++)
		{
		flag[i][j] = 0;
		assert(pic);
		switch(pic[i][j])
		{
			case 0:
			flag[i][j] = 1<<1;
			break;

			case 1:
			flag[i][j] = 1<<2;
			break;

			case 2:
			flag[i][j] = 1<<3;
			break;

			case 3:
			flag[i][j] = 1<<4;
			break;

			case 4:
			flag[i][j] = 1<<0;
			break;
		}

			if(!isfluid(pic[i][j]))
			{
				if(isfluid(pic[i+1][j]) && i<imax-1)
				{
				flag[i][j] |= 1<<8;
				}
				if(isfluid(pic[i-1][j]) && i>0)
				{
				flag[i][j] |= 1<<7;
				}
				if(isfluid(pic[i][j+1]) == 4 && j<jmax-1)
				{
				flag[i][j] |= 1<<5;
				}
				if(isfluid(pic[i][j-1]) == 4 && j>0)
				{
				flag[i][j] |= 1<<6;
				}

			}


		}

	}
	printf("Flags set using geometry file.\n");


}
