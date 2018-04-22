void calculate_rs(dt,dx,dy,imax,jmax,F,G,RS)
{

for(int i=1; i<imax; i++)
   {
	
        for(int j=1; j<jmax; j++)
		{
		
		RS[i][j] =  (1/dt)*((((F[i][j])-(F[i-1][j]))/(dx)) + (((G[i][j])-(G[i][j-1]^))/(dy)))
		
		}
		
	}

}