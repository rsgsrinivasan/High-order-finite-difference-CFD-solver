#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
	int iteration, i, j, k, gar, gar1, NUMNP;
	double u, v, w, rho, p, t, e, mu, gar2, gar3, gar4;
	FILE *fp = NULL;
	FILE *fp2 = NULL;
	char line[1000];
	iteration = 67000;
	NUMNP = 9250686;
	
	fp = fopen("nodefile_67000.plt","r");
	fp2 = fopen("restart_file.neu","w");
	while(fgets(line, 80, fp) != NULL)
	{
		i++;
		if(i >= 4)
		{
			sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &gar2, &gar3, &gar4, &u, &v, &w, &rho, &p, &t, &e, &mu, &gar, &gar1);
			fprintf(fp2,"%e %e %e %e %e %e %e %e %d\n", u, v, w, rho, p, t, e, mu, iteration);
			j++;
			if (j == NUMNP)
			{
				break;
			}
		}		
	}
	
	
	
	
}