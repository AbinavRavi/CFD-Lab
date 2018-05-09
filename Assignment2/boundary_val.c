#include <stdio.h>
#include <string.h>
#include "boundary_val.h"
#include "init.h"
#include "helper.h"

void boundaryvalues(int imax, int jmax, double **U, double **V,  int wl, int wr, int wt, int wb, int **Flag)
{
	int i, j;

	/* Iterate across both the ceiling and floor  */
	for(i = 1; i < imax+1; i++)
	{
		/* FLOOR BC*/
		if (wb == 1)	/* no-slip BC for floor*/
		{
			U[i][0] = -U[i][1];
			V[i][0]	= 0;
		}
		else if (wb == 2)	/* free-slip BC for floor*/
		{
			V[i][0]	= 0;
			U[i][0] = U[i][1];
		}
		else if (wb == 3)	/* outflow BC for floor*/
		{
			V[i][0]	= V[i][1];
			U[i][0] = U[i][1];
		}

		/* CEILING BC */
		if (wt == 1)	/* no-slip BC for ceiling*/
		{
			U[i][jmax+1] 	= -U[i][jmax];
			V[i][jmax]		= 0;
		}
		else if (wt == 2)	/* free-slip BC for ceiling*/
		{
			V[i][jmax]		= 0;
			U[i][jmax+1] 	= U[i][jmax];

		}
		else if (wt == 3) /* (wb == 3) outflow BC for ceiling*/
		{
			U[i][jmax+1] 	= U[i][jmax];
			V[i][jmax] 		= V[i][jmax-1];

		}
	}

	/* Iterate across both the walls*/
	for(j = 1; j < jmax+1; j++)
	{
		/* LEFT wall BC*/
		if (wl == 1)	/* no-slip BC for left wall*/
		{
			V[0][j]	= -V[1][j];
			U[0][j]	= 0;
		}
		else if (wl == 2)	/* free-slip BC for left wall*/
		{
			U[0][j]	= 0;
			V[0][j]	= V[1][j];
		}
		else if (wl == 3)	/* outflow BC for left wall*/
		{
			U[0][j]	= U[1][j];
			V[0][j]	= V[1][j];
		}

		/* RIGHT wall BC*/
		if (wr == 1)	/* no-slip BC for right wall*/
		{
			V[imax+1][j] 	= -V[imax][j];
			U[imax][j]		= 0;
		}
		else if (wr == 2)	/* free-slip BC for right wall*/
		{
			U[imax][j]		= 0;
			V[imax+1][j] 	= V[imax][j];
		}
		else if (wr == 3)	/* outflow BC for right wall*/
		{
			U[imax][j] 		= U[imax-1][j];
			V[imax+1][j] 	= V[imax][j];
		}
	}

	/*
		loop through the Flag array and find each cells have the flag B_N, B_W, B_O, B_S
	 */
	for(i = 1 ; i < imax + 1 ; i++)
	{
		for(j = 1 ; j < jmax + 1 ; j++)
		{
			/* check the boundary cells */
			if(Flag[i][j] == B_O)
			{
				U[i][j] = 0;
				V[i][j - 1] = - V[i + 1][j - 1];
				V[i][j]	= -V[i + 1][j];
			}
			else if(Flag[i][j] == B_N)
			{
				V[i][j] = 0;
				U[i - 1][j] = - U[i - 1][j + 1];
				U[i][j]	= -U[i][j + 1];
			}
			else if(Flag[i][j] == B_W)
			{
				U[i - 1][j] = 0;
				V[i][j - 1] = - V[i - 1][j - 1];
				V[i][j]	= -V[i - 1][j];
			}
			else if(Flag[i][j] == B_S)
			{
				V[i][j - 1] = 0;
				U[i][j] = - U[i][j - 1];
				U[i - 1][j]	= - U[i - 1][j - 1];
			}
			else if(Flag[i][j] == B_NO)
			{
				U[i][j] = 0;
				V[i][j] = 0;
				U[i - 1][j]	= -U[i - 1][j + 1];
				V[i][j - 1]	= -V[i + 1][j - 1];
			}
			else if(Flag[i][j] == B_NW)
			{
				U[i - 1][j] = 0;
				V[i][j] = 0;
				U[i][j]	= -U[i][j + 1];
				V[i][j - 1]	= -V[i - 1][j - 1];
			}
			else if(Flag[i][j] == B_SO)
			{
				U[i][j] = 0;
				V[i][j - 1] = 0;
				U[i - 1][j]	= -U[i - 1][j - 1];
				V[i][j]	= -V[i + 1][j];
			}
			else if(Flag[i][j] == B_SW)
			{
				U[i - 1][j] = 0;
				V[i][j - 1] = 0;
				U[i][j]	= -U[i][j - 1];
				V[i][j]	= -V[i - 1][j];
			}
		}
	}
}
