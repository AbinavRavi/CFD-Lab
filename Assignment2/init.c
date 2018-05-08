#include "helper.h"
#include "init.h"

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
		                double *dt_value,         /* time for output */
                                /*for Temperature inclusion/
                    double *TI,               /*Temperature*/
                    double *Pr,               /*Prandtl Number*/
                    double *beta,           /*Coefficient of Thermal Expansion*/
	            char *problem,
                    char *geometry)
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
   READ_DOUBLE( szFileName, *TI);
   READ_DOUBLE(szFileName, *Pr);
   READ_DOUBLE(szFileName, *beta);
   READ_STRING(szFileName, *problem);
   READ_STRING(szFileName, *geometry);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);
	
	void init_uvp(
              double UI,
              double VI,
              double PI,             /* from First Assignment To initialise Values*/
              int imax,
              int jmax,
              double **U,
              double **V,
              double **P
             ){

    init_matrix(U, 0, imax ,  0, jmax+1, UI);
    init_matrix(V, 0, imax+1, 0, jmax  , VI);
    init_matrix(P, 0, imax+1, 0, jmax+1, PI);

}
	
void init_flag(
  char* problem,
  char* geometry,
  int imax,
  int jmax,
  int** Flag)
{
	 int **Picture = NULL;

  Picture = read_pgm(geometry);
	
	for(int i=0; i<imax+1; i++){
    for(int j=0; j<jmax+1; j++){
      Flag[i][j] = Picture[i][j]; /* still need to be modified*/ /* No idea to Proceed */
    }
  }

}


   return 1;
}
