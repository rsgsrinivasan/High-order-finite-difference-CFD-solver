#include <stdio.h>
#include<stdlib.h> 
#include<math.h>
#include"function.h"
void writepltfile(int NELEM, int NUMNP, MNODE *node, ELEM *CD, double **u, double **v, double **w, double **rho, double **p, double **t, double **e, double **mu, int iteration, double *DUCROS)
{
int i;
FILE *fp;
char filename[100];
sprintf(filename,"nodefile_%d.dat",iteration);
fp= fopen(filename,"w");
fprintf(fp,"TITLE = \"Node file\"\n");
fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"RHO\", \"P\", \"T\", \"e\", \"mu\", \"loc\", \"Proc\", \"DUCROS\",\n");
fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK, SOLUTIONTIME=%d\n",NUMNP,NELEM,iteration);
for(i=0;i<NUMNP;i++)
	{
		fprintf(fp,"%f\t%f\t%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n", node[i+1].x, node[i+1].y, node[i+1].z, u[0][i+1], v[0][i+1], w[0][i+1], rho[0][i+1], p[0][i+1], t[0][i+1], e[0][i+1], mu[0][i+1], node[i+1].loc, node[i+1].proc, DUCROS[i+1]);
	}
	fprintf(fp,"\n\n\n");
	
	for(i=1; i<=NELEM; i++)
	{
		fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CD[i].connect[0], CD[i].connect[1], CD[i].connect[5], CD[i].connect[4], CD[i].connect[3], CD[i].connect[2], CD[i].connect[6], CD[i].connect[7] );
	}
	
	fclose(fp);
	//exit(0);
}

void restart_file(int NELEM, int NUMNP, MNODE *node, ELEM *CD, double **u, double **v, double **w, double **rho, double **p, double **t, double **e, double **mu, int iteration, int restart_num)
{
int i;
char filename[100];
sprintf(filename,"restart_file_%d.neu",restart_num);
FILE *fp = fopen(filename,"w");

for(i=0;i<NUMNP;i++)
	{
		fprintf(fp,"%e %e %e %e %e %e %e %e %d\n", u[0][i+1], v[0][i+1], w[0][i+1], rho[0][i+1], p[0][i+1], t[0][i+1], e[0][i+1], mu[0][i+1], iteration);
	}
	fprintf(fp,"\n\n\n");
		
	fclose(fp);
	//exit(0);
}

void writepltfile_jacob(int NELEM, int NUMNP, MNODE *node, ELEM *CD, double **u, double **v, double **rho, double **p, double **t, JACOB *jacobian, TRANSFORMED *metric, double *det)
{
int i;
FILE *fp = fopen("nodefilejacob.plt","w");
fprintf(fp,"TITLE = \"Node file\"\n");
fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"x_zeta\", \"x_eta\", \"y_zeta\", \"y_eta\", \"zeta_x\", \"eta_x\", \"zeta_y\", \"eta_y\", \"det\", \n");
fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",NUMNP,NELEM);
for(i=0;i<NUMNP;i++)
	{
		fprintf(fp,"%f\t%f\t%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", node[i+1].x, node[i+1].y, node[i+1].z, jacobian[i+1].x_zeta, jacobian[i+1].x_eta, jacobian[i+1].y_zeta, jacobian[i+1].y_eta, metric[i+1].zeta_x, metric[i+1].eta_x, metric[i+1].zeta_y, metric[i+1].eta_y, det[i+1] );
	}
	fprintf(fp,"\n\n\n");
	
	for(i=1; i<=NELEM; i++)
	{
		fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CD[i].connect[0], CD[i].connect[1], CD[i].connect[5], CD[i].connect[4], CD[i].connect[3], CD[i].connect[2], CD[i].connect[6], CD[i].connect[7] );
	}
	
	fclose(fp);
	//exit(0);
}





