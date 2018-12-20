/************************************************************************************************************************************************************************
								Lax-Friedrich
								Flux splitting
**************************************************************************************************************************************************************************/

#include <stdio.h>
#include<stdlib.h> 
#include<math.h>
#include"function.h"

void weno_solver_0(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus)
{				
	int k,m;	
	double det_a_f;
	double eigen_E[5];
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi, eta_xi, eta_yi, eta_zi, xi_xi, xi_yi, xi_zi;
	
	roe_a[i] = sqrt(0.4*(roe_h[i]-0.5*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i])));	

//	if (hk == 0)
	{
		zeta_xi = metric[i].eta_xjp;
		zeta_yi = metric[i].eta_yjp;
		zeta_zi = metric[i].eta_zjp;
	}
	
	kk = sqrt(zeta_xi*zeta_xi+zeta_yi*zeta_yi+zeta_zi*zeta_zi);
	
	kx_bar = zeta_xi/kk;
	ky_bar = zeta_yi/kk;
	kz_bar = zeta_zi/kk;

	theta = kx_bar*roe_u[i]+ky_bar*roe_v[i]+kz_bar*roe_w[i];
	phi_sq = 0.5*0.4*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i]);
	alpha = roe_rho[i]/(sqrt(2.0)*roe_a[i]);
	beta = 1.0/(sqrt(2.0)*roe_rho[i]*roe_a[i]);
	
	r_eigen_f[0][0]= kx_bar;
	r_eigen_f[0][1]= ky_bar;
	r_eigen_f[0][2]= kz_bar;
	r_eigen_f[0][3]= alpha;
	r_eigen_f[0][4]= alpha;

	r_eigen_f[1][0]= kx_bar*roe_u[i];
	r_eigen_f[1][1]= ky_bar*roe_u[i]-kz_bar*roe_rho[i];
	r_eigen_f[1][2]= kz_bar*roe_u[i]+ky_bar*roe_rho[i];
	r_eigen_f[1][3]= alpha*(roe_u[i]+kx_bar*roe_a[i]);
	r_eigen_f[1][4]= alpha*(roe_u[i]-kx_bar*roe_a[i]);

	r_eigen_f[2][0]= kx_bar*roe_v[i]+kz_bar*roe_rho[i];
	r_eigen_f[2][1]= ky_bar*roe_v[i];
	r_eigen_f[2][2]= kz_bar*roe_v[i]-kx_bar*roe_rho[i];
	r_eigen_f[2][3]= alpha*(roe_v[i]+ky_bar*roe_a[i]);
	r_eigen_f[2][4]= alpha*(roe_v[i]-ky_bar*roe_a[i]);
	
	r_eigen_f[3][0]= kx_bar*roe_w[i]-ky_bar*roe_rho[i];
	r_eigen_f[3][1]= ky_bar*roe_w[i]+kx_bar*roe_rho[i];
	r_eigen_f[3][2]= kz_bar*roe_w[i];
	r_eigen_f[3][3]= alpha*(roe_w[i]+kz_bar*roe_a[i]);
	r_eigen_f[3][4]= alpha*(roe_w[i]-kz_bar*roe_a[i]);
	
	r_eigen_f[4][0]= ((kx_bar*phi_sq)/0.4)+roe_rho[i]*(kz_bar*roe_v[i]-ky_bar*roe_w[i]);
	r_eigen_f[4][1]= ((ky_bar*phi_sq)/0.4)+roe_rho[i]*(kx_bar*roe_w[i]-kz_bar*roe_u[i]);
	r_eigen_f[4][2]= ((kz_bar*phi_sq)/0.4)+roe_rho[i]*(ky_bar*roe_u[i]-kx_bar*roe_v[i]);
	r_eigen_f[4][3]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))+theta*roe_a[i]);
	r_eigen_f[4][4]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))-theta*roe_a[i]);
	
	/*********************************************/
	
	l_eigen_f[0][0]= kx_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kz_bar*roe_v[i]-ky_bar*roe_w[i])/roe_rho[i]);
	l_eigen_f[0][1]= kx_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[0][2]= (kx_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]))+(kz_bar/roe_rho[i]);
	l_eigen_f[0][3]= (kx_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))-(ky_bar/roe_rho[i]);
	l_eigen_f[0][4]= -kx_bar*0.4/(roe_a[i]*roe_a[i]);
			
	l_eigen_f[1][0]= ky_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kx_bar*roe_w[i]-kz_bar*roe_u[i])/roe_rho[i]);
	l_eigen_f[1][1]= (ky_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]))-kz_bar/roe_rho[i];
	l_eigen_f[1][2]= ky_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[1][3]= (ky_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))+(kx_bar/roe_rho[i]);
	l_eigen_f[1][4]=	-ky_bar*0.4/(roe_a[i]*roe_a[i]);
				
	l_eigen_f[2][0]= kz_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((ky_bar*roe_u[i]-kx_bar*roe_v[i])/roe_rho[i]);
	l_eigen_f[2][1]= kz_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i])+(ky_bar/roe_rho[i]);
	l_eigen_f[2][2]= kz_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i])-(kx_bar/roe_rho[i]);
	l_eigen_f[2][3]= kz_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[2][4]= -kz_bar*0.4/(roe_a[i]*roe_a[i]);
	
	l_eigen_f[3][0]= beta*(phi_sq-theta*roe_a[i]);
	l_eigen_f[3][1]= -beta*(0.4*roe_u[i]-kx_bar*roe_a[i]);
	l_eigen_f[3][2]= -beta*(0.4*roe_v[i]-ky_bar*roe_a[i]);
	l_eigen_f[3][3]= -beta*(0.4*roe_w[i]-kz_bar*roe_a[i]);
	l_eigen_f[3][4]= beta*0.4;

	l_eigen_f[4][0]= beta*(phi_sq+theta*roe_a[i]);
	l_eigen_f[4][1]= -beta*(0.4*roe_u[i]+kx_bar*roe_a[i]);
	l_eigen_f[4][2]= -beta*(0.4*roe_v[i]+ky_bar*roe_a[i]);
	l_eigen_f[4][3]= -beta*(0.4*roe_w[i]+kz_bar*roe_a[i]);
	l_eigen_f[4][4]= beta*0.4;
	
	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/			

	m=0;
	for(k=0; k<5; k++)
	{	
	//	if (hk == 0)
		{
			f_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][4].ip;			
			f_iplus[node[node[i].n_n[0]].n_n[0]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[0]].n_n[0]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[0]].n_n[0]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[0]].n_n[0]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[0]].n_n[0]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[0]].n_n[0]][4].ip;			
			f_iplus[node[i].n_n[0]][k] = l_eigen_f[m][0]*Q[node[i].n_n[0]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[0]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[0]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[0]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[0]][4].ip;			
			f_iplus[i][k] = l_eigen_f[m][0]*Q[i][0].ip+l_eigen_f[m][1]*Q[i][1].ip+l_eigen_f[m][2]*Q[i][2].ip+l_eigen_f[m][3]*Q[i][3].ip+l_eigen_f[m][4]*Q[i][4].ip;
			f_iplus[node[i].n_n[2]][k] = l_eigen_f[m][0]*Q[node[i].n_n[2]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[2]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[2]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[2]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[2]][4].ip;
			f_iplus[node[node[i].n_n[2]].n_n[2]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[2]].n_n[2]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[2]].n_n[2]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[2]].n_n[2]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[2]].n_n[2]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[2]].n_n[2]][4].ip;			
			f_iplus[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][4].ip;						
		}
		
		m++; 
	}
}


void weno_solver_1(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus)
{				
	int k,m;	
	double det_a_f;
	double eigen_E[5];
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi, eta_xi, eta_yi, eta_zi, xi_xi, xi_yi, xi_zi;
	
	roe_a[i] = sqrt(0.4*(roe_h[i]-0.5*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i])));	

	zeta_xi = metric[i].zeta_xip;
	zeta_yi = metric[i].zeta_yip;
	zeta_zi = metric[i].zeta_zip;
	
	kk = sqrt(zeta_xi*zeta_xi+zeta_yi*zeta_yi+zeta_zi*zeta_zi);
	
	kx_bar = zeta_xi/kk;
	ky_bar = zeta_yi/kk;
	kz_bar = zeta_zi/kk;

	theta = kx_bar*roe_u[i]+ky_bar*roe_v[i]+kz_bar*roe_w[i];
	phi_sq = 0.5*0.4*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i]);
	alpha = roe_rho[i]/(sqrt(2.0)*roe_a[i]);
	beta = 1.0/(sqrt(2.0)*roe_rho[i]*roe_a[i]);
	
	r_eigen_f[0][0]= kx_bar;
	r_eigen_f[0][1]= ky_bar;
	r_eigen_f[0][2]= kz_bar;
	r_eigen_f[0][3]= alpha;
	r_eigen_f[0][4]= alpha;

	r_eigen_f[1][0]= kx_bar*roe_u[i];
	r_eigen_f[1][1]= ky_bar*roe_u[i]-kz_bar*roe_rho[i];
	r_eigen_f[1][2]= kz_bar*roe_u[i]+ky_bar*roe_rho[i];
	r_eigen_f[1][3]= alpha*(roe_u[i]+kx_bar*roe_a[i]);
	r_eigen_f[1][4]= alpha*(roe_u[i]-kx_bar*roe_a[i]);

	r_eigen_f[2][0]= kx_bar*roe_v[i]+kz_bar*roe_rho[i];
	r_eigen_f[2][1]= ky_bar*roe_v[i];
	r_eigen_f[2][2]= kz_bar*roe_v[i]-kx_bar*roe_rho[i];
	r_eigen_f[2][3]= alpha*(roe_v[i]+ky_bar*roe_a[i]);
	r_eigen_f[2][4]= alpha*(roe_v[i]-ky_bar*roe_a[i]);
	
	r_eigen_f[3][0]= kx_bar*roe_w[i]-ky_bar*roe_rho[i];
	r_eigen_f[3][1]= ky_bar*roe_w[i]+kx_bar*roe_rho[i];
	r_eigen_f[3][2]= kz_bar*roe_w[i];
	r_eigen_f[3][3]= alpha*(roe_w[i]+kz_bar*roe_a[i]);
	r_eigen_f[3][4]= alpha*(roe_w[i]-kz_bar*roe_a[i]);
	
	r_eigen_f[4][0]= ((kx_bar*phi_sq)/0.4)+roe_rho[i]*(kz_bar*roe_v[i]-ky_bar*roe_w[i]);
	r_eigen_f[4][1]= ((ky_bar*phi_sq)/0.4)+roe_rho[i]*(kx_bar*roe_w[i]-kz_bar*roe_u[i]);
	r_eigen_f[4][2]= ((kz_bar*phi_sq)/0.4)+roe_rho[i]*(ky_bar*roe_u[i]-kx_bar*roe_v[i]);
	r_eigen_f[4][3]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))+theta*roe_a[i]);
	r_eigen_f[4][4]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))-theta*roe_a[i]);
	
	/*********************************************/
	
	l_eigen_f[0][0]= kx_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kz_bar*roe_v[i]-ky_bar*roe_w[i])/roe_rho[i]);
	l_eigen_f[0][1]= kx_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[0][2]= (kx_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]))+(kz_bar/roe_rho[i]);
	l_eigen_f[0][3]= (kx_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))-(ky_bar/roe_rho[i]);
	l_eigen_f[0][4]= -kx_bar*0.4/(roe_a[i]*roe_a[i]);
			
	l_eigen_f[1][0]= ky_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kx_bar*roe_w[i]-kz_bar*roe_u[i])/roe_rho[i]);
	l_eigen_f[1][1]= (ky_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]))-kz_bar/roe_rho[i];
	l_eigen_f[1][2]= ky_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[1][3]= (ky_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))+(kx_bar/roe_rho[i]);
	l_eigen_f[1][4]=	-ky_bar*0.4/(roe_a[i]*roe_a[i]);
				
	l_eigen_f[2][0]= kz_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((ky_bar*roe_u[i]-kx_bar*roe_v[i])/roe_rho[i]);
	l_eigen_f[2][1]= kz_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i])+(ky_bar/roe_rho[i]);
	l_eigen_f[2][2]= kz_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i])-(kx_bar/roe_rho[i]);
	l_eigen_f[2][3]= kz_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[2][4]= -kz_bar*0.4/(roe_a[i]*roe_a[i]);
	
	l_eigen_f[3][0]= beta*(phi_sq-theta*roe_a[i]);
	l_eigen_f[3][1]= -beta*(0.4*roe_u[i]-kx_bar*roe_a[i]);
	l_eigen_f[3][2]= -beta*(0.4*roe_v[i]-ky_bar*roe_a[i]);
	l_eigen_f[3][3]= -beta*(0.4*roe_w[i]-kz_bar*roe_a[i]);
	l_eigen_f[3][4]= beta*0.4;

	l_eigen_f[4][0]= beta*(phi_sq+theta*roe_a[i]);
	l_eigen_f[4][1]= -beta*(0.4*roe_u[i]+kx_bar*roe_a[i]);
	l_eigen_f[4][2]= -beta*(0.4*roe_v[i]+ky_bar*roe_a[i]);
	l_eigen_f[4][3]= -beta*(0.4*roe_w[i]+kz_bar*roe_a[i]);
	l_eigen_f[4][4]= beta*0.4;
	
	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/			

	m=0;
	for(k=0; k<5; k++)
	{	
		f_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][4].ip;			
		f_iplus[node[node[i].n_n[1]].n_n[1]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[1]].n_n[1]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[1]].n_n[1]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[1]].n_n[1]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[1]].n_n[1]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[1]].n_n[1]][4].ip;			
		f_iplus[node[i].n_n[1]][k] = l_eigen_f[m][0]*Q[node[i].n_n[1]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[1]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[1]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[1]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[1]][4].ip;			
		f_iplus[i][k] = l_eigen_f[m][0]*Q[i][0].ip+l_eigen_f[m][1]*Q[i][1].ip+l_eigen_f[m][2]*Q[i][2].ip+l_eigen_f[m][3]*Q[i][3].ip+l_eigen_f[m][4]*Q[i][4].ip;
		f_iplus[node[i].n_n[3]][k] = l_eigen_f[m][0]*Q[node[i].n_n[3]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[3]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[3]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[3]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[3]][4].ip;
		f_iplus[node[node[i].n_n[3]].n_n[3]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[3]].n_n[3]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[3]].n_n[3]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[3]].n_n[3]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[3]].n_n[3]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[3]].n_n[3]][4].ip;			
		f_iplus[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][4].ip;						

		
		m++; 
	}
}

void weno_solver_2(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus)
{				
	int k,m;	
	double det_a_f;
	double eigen_E[5];
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi, eta_xi, eta_yi, eta_zi, xi_xi, xi_yi, xi_zi;
	
	roe_a[i] = sqrt(0.4*(roe_h[i]-0.5*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i])));	

	zeta_xi = metric[i].eta_xjm;
	zeta_yi = metric[i].eta_yjm;
	zeta_zi = metric[i].eta_zjm;	
	
	
	kk = sqrt(zeta_xi*zeta_xi+zeta_yi*zeta_yi+zeta_zi*zeta_zi);
	
	kx_bar = zeta_xi/kk;
	ky_bar = zeta_yi/kk;
	kz_bar = zeta_zi/kk;

	theta = kx_bar*roe_u[i]+ky_bar*roe_v[i]+kz_bar*roe_w[i];
	phi_sq = 0.5*0.4*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i]);
	alpha = roe_rho[i]/(sqrt(2.0)*roe_a[i]);
	beta = 1.0/(sqrt(2.0)*roe_rho[i]*roe_a[i]);
	
	r_eigen_f[0][0]= kx_bar;
	r_eigen_f[0][1]= ky_bar;
	r_eigen_f[0][2]= kz_bar;
	r_eigen_f[0][3]= alpha;
	r_eigen_f[0][4]= alpha;

	r_eigen_f[1][0]= kx_bar*roe_u[i];
	r_eigen_f[1][1]= ky_bar*roe_u[i]-kz_bar*roe_rho[i];
	r_eigen_f[1][2]= kz_bar*roe_u[i]+ky_bar*roe_rho[i];
	r_eigen_f[1][3]= alpha*(roe_u[i]+kx_bar*roe_a[i]);
	r_eigen_f[1][4]= alpha*(roe_u[i]-kx_bar*roe_a[i]);

	r_eigen_f[2][0]= kx_bar*roe_v[i]+kz_bar*roe_rho[i];
	r_eigen_f[2][1]= ky_bar*roe_v[i];
	r_eigen_f[2][2]= kz_bar*roe_v[i]-kx_bar*roe_rho[i];
	r_eigen_f[2][3]= alpha*(roe_v[i]+ky_bar*roe_a[i]);
	r_eigen_f[2][4]= alpha*(roe_v[i]-ky_bar*roe_a[i]);
	
	r_eigen_f[3][0]= kx_bar*roe_w[i]-ky_bar*roe_rho[i];
	r_eigen_f[3][1]= ky_bar*roe_w[i]+kx_bar*roe_rho[i];
	r_eigen_f[3][2]= kz_bar*roe_w[i];
	r_eigen_f[3][3]= alpha*(roe_w[i]+kz_bar*roe_a[i]);
	r_eigen_f[3][4]= alpha*(roe_w[i]-kz_bar*roe_a[i]);
	
	r_eigen_f[4][0]= ((kx_bar*phi_sq)/0.4)+roe_rho[i]*(kz_bar*roe_v[i]-ky_bar*roe_w[i]);
	r_eigen_f[4][1]= ((ky_bar*phi_sq)/0.4)+roe_rho[i]*(kx_bar*roe_w[i]-kz_bar*roe_u[i]);
	r_eigen_f[4][2]= ((kz_bar*phi_sq)/0.4)+roe_rho[i]*(ky_bar*roe_u[i]-kx_bar*roe_v[i]);
	r_eigen_f[4][3]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))+theta*roe_a[i]);
	r_eigen_f[4][4]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))-theta*roe_a[i]);
	
	/*********************************************/
	
	l_eigen_f[0][0]= kx_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kz_bar*roe_v[i]-ky_bar*roe_w[i])/roe_rho[i]);
	l_eigen_f[0][1]= kx_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[0][2]= (kx_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]))+(kz_bar/roe_rho[i]);
	l_eigen_f[0][3]= (kx_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))-(ky_bar/roe_rho[i]);
	l_eigen_f[0][4]= -kx_bar*0.4/(roe_a[i]*roe_a[i]);
			
	l_eigen_f[1][0]= ky_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kx_bar*roe_w[i]-kz_bar*roe_u[i])/roe_rho[i]);
	l_eigen_f[1][1]= (ky_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]))-kz_bar/roe_rho[i];
	l_eigen_f[1][2]= ky_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[1][3]= (ky_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))+(kx_bar/roe_rho[i]);
	l_eigen_f[1][4]=	-ky_bar*0.4/(roe_a[i]*roe_a[i]);
				
	l_eigen_f[2][0]= kz_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((ky_bar*roe_u[i]-kx_bar*roe_v[i])/roe_rho[i]);
	l_eigen_f[2][1]= kz_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i])+(ky_bar/roe_rho[i]);
	l_eigen_f[2][2]= kz_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i])-(kx_bar/roe_rho[i]);
	l_eigen_f[2][3]= kz_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[2][4]= -kz_bar*0.4/(roe_a[i]*roe_a[i]);
	
	l_eigen_f[3][0]= beta*(phi_sq-theta*roe_a[i]);
	l_eigen_f[3][1]= -beta*(0.4*roe_u[i]-kx_bar*roe_a[i]);
	l_eigen_f[3][2]= -beta*(0.4*roe_v[i]-ky_bar*roe_a[i]);
	l_eigen_f[3][3]= -beta*(0.4*roe_w[i]-kz_bar*roe_a[i]);
	l_eigen_f[3][4]= beta*0.4;

	l_eigen_f[4][0]= beta*(phi_sq+theta*roe_a[i]);
	l_eigen_f[4][1]= -beta*(0.4*roe_u[i]+kx_bar*roe_a[i]);
	l_eigen_f[4][2]= -beta*(0.4*roe_v[i]+ky_bar*roe_a[i]);
	l_eigen_f[4][3]= -beta*(0.4*roe_w[i]+kz_bar*roe_a[i]);
	l_eigen_f[4][4]= beta*0.4;
	
	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/			

	m=0;
	for(k=0; k<5; k++)
	{	
	//	if (hk == 2)
		{
			f_iminus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][4].ip;			
			f_iminus[node[node[i].n_n[0]].n_n[0]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[0]].n_n[0]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[0]].n_n[0]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[0]].n_n[0]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[0]].n_n[0]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[0]].n_n[0]][4].ip;			
			f_iminus[node[i].n_n[0]][k] = l_eigen_f[m][0]*Q[node[i].n_n[0]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[0]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[0]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[0]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[0]][4].ip;			
			f_iminus[i][k] = l_eigen_f[m][0]*Q[i][0].ip+l_eigen_f[m][1]*Q[i][1].ip+l_eigen_f[m][2]*Q[i][2].ip+l_eigen_f[m][3]*Q[i][3].ip+l_eigen_f[m][4]*Q[i][4].ip;
			f_iminus[node[i].n_n[2]][k] = l_eigen_f[m][0]*Q[node[i].n_n[2]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[2]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[2]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[2]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[2]][4].ip;
			f_iminus[node[node[i].n_n[2]].n_n[2]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[2]].n_n[2]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[2]].n_n[2]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[2]].n_n[2]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[2]].n_n[2]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[2]].n_n[2]][4].ip;			
			f_iminus[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][4].ip;						
		}
		
		m++; 
	}
}

void weno_solver_3(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus)
{				
	int k,m;	
	double det_a_f;
	double eigen_E[5];
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi, eta_xi, eta_yi, eta_zi, xi_xi, xi_yi, xi_zi;
	
	roe_a[i] = sqrt(0.4*(roe_h[i]-0.5*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i])));	

//	if(hk == 3)
	{
		zeta_xi = metric[i].zeta_xim;
		zeta_yi = metric[i].zeta_yim;
		zeta_zi = metric[i].zeta_zim;
	}
	
	kk = sqrt(zeta_xi*zeta_xi+zeta_yi*zeta_yi+zeta_zi*zeta_zi);
	
	kx_bar = zeta_xi/kk;
	ky_bar = zeta_yi/kk;
	kz_bar = zeta_zi/kk;

	theta = kx_bar*roe_u[i]+ky_bar*roe_v[i]+kz_bar*roe_w[i];
	phi_sq = 0.5*0.4*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i]);
	alpha = roe_rho[i]/(sqrt(2.0)*roe_a[i]);
	beta = 1.0/(sqrt(2.0)*roe_rho[i]*roe_a[i]);
	
	r_eigen_f[0][0]= kx_bar;
	r_eigen_f[0][1]= ky_bar;
	r_eigen_f[0][2]= kz_bar;
	r_eigen_f[0][3]= alpha;
	r_eigen_f[0][4]= alpha;

	r_eigen_f[1][0]= kx_bar*roe_u[i];
	r_eigen_f[1][1]= ky_bar*roe_u[i]-kz_bar*roe_rho[i];
	r_eigen_f[1][2]= kz_bar*roe_u[i]+ky_bar*roe_rho[i];
	r_eigen_f[1][3]= alpha*(roe_u[i]+kx_bar*roe_a[i]);
	r_eigen_f[1][4]= alpha*(roe_u[i]-kx_bar*roe_a[i]);

	r_eigen_f[2][0]= kx_bar*roe_v[i]+kz_bar*roe_rho[i];
	r_eigen_f[2][1]= ky_bar*roe_v[i];
	r_eigen_f[2][2]= kz_bar*roe_v[i]-kx_bar*roe_rho[i];
	r_eigen_f[2][3]= alpha*(roe_v[i]+ky_bar*roe_a[i]);
	r_eigen_f[2][4]= alpha*(roe_v[i]-ky_bar*roe_a[i]);
	
	r_eigen_f[3][0]= kx_bar*roe_w[i]-ky_bar*roe_rho[i];
	r_eigen_f[3][1]= ky_bar*roe_w[i]+kx_bar*roe_rho[i];
	r_eigen_f[3][2]= kz_bar*roe_w[i];
	r_eigen_f[3][3]= alpha*(roe_w[i]+kz_bar*roe_a[i]);
	r_eigen_f[3][4]= alpha*(roe_w[i]-kz_bar*roe_a[i]);
	
	r_eigen_f[4][0]= ((kx_bar*phi_sq)/0.4)+roe_rho[i]*(kz_bar*roe_v[i]-ky_bar*roe_w[i]);
	r_eigen_f[4][1]= ((ky_bar*phi_sq)/0.4)+roe_rho[i]*(kx_bar*roe_w[i]-kz_bar*roe_u[i]);
	r_eigen_f[4][2]= ((kz_bar*phi_sq)/0.4)+roe_rho[i]*(ky_bar*roe_u[i]-kx_bar*roe_v[i]);
	r_eigen_f[4][3]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))+theta*roe_a[i]);
	r_eigen_f[4][4]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))-theta*roe_a[i]);
	
	/*********************************************/
	
	l_eigen_f[0][0]= kx_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kz_bar*roe_v[i]-ky_bar*roe_w[i])/roe_rho[i]);
	l_eigen_f[0][1]= kx_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[0][2]= (kx_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]))+(kz_bar/roe_rho[i]);
	l_eigen_f[0][3]= (kx_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))-(ky_bar/roe_rho[i]);
	l_eigen_f[0][4]= -kx_bar*0.4/(roe_a[i]*roe_a[i]);
			
	l_eigen_f[1][0]= ky_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kx_bar*roe_w[i]-kz_bar*roe_u[i])/roe_rho[i]);
	l_eigen_f[1][1]= (ky_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]))-kz_bar/roe_rho[i];
	l_eigen_f[1][2]= ky_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[1][3]= (ky_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))+(kx_bar/roe_rho[i]);
	l_eigen_f[1][4]=	-ky_bar*0.4/(roe_a[i]*roe_a[i]);
				
	l_eigen_f[2][0]= kz_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((ky_bar*roe_u[i]-kx_bar*roe_v[i])/roe_rho[i]);
	l_eigen_f[2][1]= kz_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i])+(ky_bar/roe_rho[i]);
	l_eigen_f[2][2]= kz_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i])-(kx_bar/roe_rho[i]);
	l_eigen_f[2][3]= kz_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[2][4]= -kz_bar*0.4/(roe_a[i]*roe_a[i]);
	
	l_eigen_f[3][0]= beta*(phi_sq-theta*roe_a[i]);
	l_eigen_f[3][1]= -beta*(0.4*roe_u[i]-kx_bar*roe_a[i]);
	l_eigen_f[3][2]= -beta*(0.4*roe_v[i]-ky_bar*roe_a[i]);
	l_eigen_f[3][3]= -beta*(0.4*roe_w[i]-kz_bar*roe_a[i]);
	l_eigen_f[3][4]= beta*0.4;

	l_eigen_f[4][0]= beta*(phi_sq+theta*roe_a[i]);
	l_eigen_f[4][1]= -beta*(0.4*roe_u[i]+kx_bar*roe_a[i]);
	l_eigen_f[4][2]= -beta*(0.4*roe_v[i]+ky_bar*roe_a[i]);
	l_eigen_f[4][3]= -beta*(0.4*roe_w[i]+kz_bar*roe_a[i]);
	l_eigen_f[4][4]= beta*0.4;
	
	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/			

	m=0;
	for(k=0; k<5; k++)
	{	
	//	if(hk == 3 )
		{
			f_iminus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][4].ip;			
			f_iminus[node[node[i].n_n[1]].n_n[1]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[1]].n_n[1]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[1]].n_n[1]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[1]].n_n[1]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[1]].n_n[1]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[1]].n_n[1]][4].ip;			
			f_iminus[node[i].n_n[1]][k] = l_eigen_f[m][0]*Q[node[i].n_n[1]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[1]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[1]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[1]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[1]][4].ip;			
			f_iminus[i][k] = l_eigen_f[m][0]*Q[i][0].ip+l_eigen_f[m][1]*Q[i][1].ip+l_eigen_f[m][2]*Q[i][2].ip+l_eigen_f[m][3]*Q[i][3].ip+l_eigen_f[m][4]*Q[i][4].ip;
			f_iminus[node[i].n_n[3]][k] = l_eigen_f[m][0]*Q[node[i].n_n[3]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[3]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[3]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[3]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[3]][4].ip;
			f_iminus[node[node[i].n_n[3]].n_n[3]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[3]].n_n[3]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[3]].n_n[3]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[3]].n_n[3]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[3]].n_n[3]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[3]].n_n[3]][4].ip;			
			f_iminus[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][4].ip;						
		}		
		m++; 
	}
}

void weno_solver_4(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus)
{				
	int k,m;	
	double det_a_f;
	double eigen_E[5];
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi, eta_xi, eta_yi, eta_zi, xi_xi, xi_yi, xi_zi;
	
	roe_a[i] = sqrt(0.4*(roe_h[i]-0.5*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i])));	

//	if (hk == 4)
	{
		zeta_xi = metric[i].xi_xkp;
		zeta_yi = metric[i].xi_ykp;
		zeta_zi = metric[i].xi_zkp;
	}
	
	kk = sqrt(zeta_xi*zeta_xi+zeta_yi*zeta_yi+zeta_zi*zeta_zi);
	
	kx_bar = zeta_xi/kk;
	ky_bar = zeta_yi/kk;
	kz_bar = zeta_zi/kk;

	theta = kx_bar*roe_u[i]+ky_bar*roe_v[i]+kz_bar*roe_w[i];
	phi_sq = 0.5*0.4*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i]);
	alpha = roe_rho[i]/(sqrt(2.0)*roe_a[i]);
	beta = 1.0/(sqrt(2.0)*roe_rho[i]*roe_a[i]);
	
	r_eigen_f[0][0]= kx_bar;
	r_eigen_f[0][1]= ky_bar;
	r_eigen_f[0][2]= kz_bar;
	r_eigen_f[0][3]= alpha;
	r_eigen_f[0][4]= alpha;

	r_eigen_f[1][0]= kx_bar*roe_u[i];
	r_eigen_f[1][1]= ky_bar*roe_u[i]-kz_bar*roe_rho[i];
	r_eigen_f[1][2]= kz_bar*roe_u[i]+ky_bar*roe_rho[i];
	r_eigen_f[1][3]= alpha*(roe_u[i]+kx_bar*roe_a[i]);
	r_eigen_f[1][4]= alpha*(roe_u[i]-kx_bar*roe_a[i]);

	r_eigen_f[2][0]= kx_bar*roe_v[i]+kz_bar*roe_rho[i];
	r_eigen_f[2][1]= ky_bar*roe_v[i];
	r_eigen_f[2][2]= kz_bar*roe_v[i]-kx_bar*roe_rho[i];
	r_eigen_f[2][3]= alpha*(roe_v[i]+ky_bar*roe_a[i]);
	r_eigen_f[2][4]= alpha*(roe_v[i]-ky_bar*roe_a[i]);
	
	r_eigen_f[3][0]= kx_bar*roe_w[i]-ky_bar*roe_rho[i];
	r_eigen_f[3][1]= ky_bar*roe_w[i]+kx_bar*roe_rho[i];
	r_eigen_f[3][2]= kz_bar*roe_w[i];
	r_eigen_f[3][3]= alpha*(roe_w[i]+kz_bar*roe_a[i]);
	r_eigen_f[3][4]= alpha*(roe_w[i]-kz_bar*roe_a[i]);
	
	r_eigen_f[4][0]= ((kx_bar*phi_sq)/0.4)+roe_rho[i]*(kz_bar*roe_v[i]-ky_bar*roe_w[i]);
	r_eigen_f[4][1]= ((ky_bar*phi_sq)/0.4)+roe_rho[i]*(kx_bar*roe_w[i]-kz_bar*roe_u[i]);
	r_eigen_f[4][2]= ((kz_bar*phi_sq)/0.4)+roe_rho[i]*(ky_bar*roe_u[i]-kx_bar*roe_v[i]);
	r_eigen_f[4][3]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))+theta*roe_a[i]);
	r_eigen_f[4][4]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))-theta*roe_a[i]);
	
	/*********************************************/
	
	l_eigen_f[0][0]= kx_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kz_bar*roe_v[i]-ky_bar*roe_w[i])/roe_rho[i]);
	l_eigen_f[0][1]= kx_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[0][2]= (kx_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]))+(kz_bar/roe_rho[i]);
	l_eigen_f[0][3]= (kx_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))-(ky_bar/roe_rho[i]);
	l_eigen_f[0][4]= -kx_bar*0.4/(roe_a[i]*roe_a[i]);
			
	l_eigen_f[1][0]= ky_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kx_bar*roe_w[i]-kz_bar*roe_u[i])/roe_rho[i]);
	l_eigen_f[1][1]= (ky_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]))-kz_bar/roe_rho[i];
	l_eigen_f[1][2]= ky_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[1][3]= (ky_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))+(kx_bar/roe_rho[i]);
	l_eigen_f[1][4]=	-ky_bar*0.4/(roe_a[i]*roe_a[i]);
				
	l_eigen_f[2][0]= kz_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((ky_bar*roe_u[i]-kx_bar*roe_v[i])/roe_rho[i]);
	l_eigen_f[2][1]= kz_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i])+(ky_bar/roe_rho[i]);
	l_eigen_f[2][2]= kz_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i])-(kx_bar/roe_rho[i]);
	l_eigen_f[2][3]= kz_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[2][4]= -kz_bar*0.4/(roe_a[i]*roe_a[i]);
	
	l_eigen_f[3][0]= beta*(phi_sq-theta*roe_a[i]);
	l_eigen_f[3][1]= -beta*(0.4*roe_u[i]-kx_bar*roe_a[i]);
	l_eigen_f[3][2]= -beta*(0.4*roe_v[i]-ky_bar*roe_a[i]);
	l_eigen_f[3][3]= -beta*(0.4*roe_w[i]-kz_bar*roe_a[i]);
	l_eigen_f[3][4]= beta*0.4;

	l_eigen_f[4][0]= beta*(phi_sq+theta*roe_a[i]);
	l_eigen_f[4][1]= -beta*(0.4*roe_u[i]+kx_bar*roe_a[i]);
	l_eigen_f[4][2]= -beta*(0.4*roe_v[i]+ky_bar*roe_a[i]);
	l_eigen_f[4][3]= -beta*(0.4*roe_w[i]+kz_bar*roe_a[i]);
	l_eigen_f[4][4]= beta*0.4;
	
	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/			

	m=0;
	for(k=0; k<5; k++)
	{	
	//	if (hk == 4)
		{
			f_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][4].ip;			
			f_iplus[node[node[i].n_n[4]].n_n[4]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[4]].n_n[4]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[4]].n_n[4]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[4]].n_n[4]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[4]].n_n[4]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[4]].n_n[4]][4].ip;			
			f_iplus[node[i].n_n[4]][k] = l_eigen_f[m][0]*Q[node[i].n_n[4]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[4]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[4]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[4]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[4]][4].ip;			
			f_iplus[i][k] = l_eigen_f[m][0]*Q[i][0].ip+l_eigen_f[m][1]*Q[i][1].ip+l_eigen_f[m][2]*Q[i][2].ip+l_eigen_f[m][3]*Q[i][3].ip+l_eigen_f[m][4]*Q[i][4].ip;
			f_iplus[node[i].n_n[5]][k] = l_eigen_f[m][0]*Q[node[i].n_n[5]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[5]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[5]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[5]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[5]][4].ip;
			f_iplus[node[node[i].n_n[5]].n_n[5]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[5]].n_n[5]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[5]].n_n[5]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[5]].n_n[5]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[5]].n_n[5]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[5]].n_n[5]][4].ip;			
			f_iplus[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][4].ip;						
		}
		
		m++; 
	}
}


void weno_solver_5(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus)
{				
	int k,m;	
	double det_a_f;
	double eigen_E[5];
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi, eta_xi, eta_yi, eta_zi, xi_xi, xi_yi, xi_zi;
	
	roe_a[i] = sqrt(0.4*(roe_h[i]-0.5*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i])));	

//	if(hk == 5)
	{
		zeta_xi = metric[i].xi_xkm;
		zeta_yi = metric[i].xi_ykm;
		zeta_zi = metric[i].xi_zkm;	
	}
	
	kk = sqrt(zeta_xi*zeta_xi+zeta_yi*zeta_yi+zeta_zi*zeta_zi);
	
	kx_bar = zeta_xi/kk;
	ky_bar = zeta_yi/kk;
	kz_bar = zeta_zi/kk;

	theta = kx_bar*roe_u[i]+ky_bar*roe_v[i]+kz_bar*roe_w[i];
	phi_sq = 0.5*0.4*(roe_u[i]*roe_u[i]+roe_v[i]*roe_v[i]+roe_w[i]*roe_w[i]);
	alpha = roe_rho[i]/(sqrt(2.0)*roe_a[i]);
	beta = 1.0/(sqrt(2.0)*roe_rho[i]*roe_a[i]);
	
	r_eigen_f[0][0]= kx_bar;
	r_eigen_f[0][1]= ky_bar;
	r_eigen_f[0][2]= kz_bar;
	r_eigen_f[0][3]= alpha;
	r_eigen_f[0][4]= alpha;

	r_eigen_f[1][0]= kx_bar*roe_u[i];
	r_eigen_f[1][1]= ky_bar*roe_u[i]-kz_bar*roe_rho[i];
	r_eigen_f[1][2]= kz_bar*roe_u[i]+ky_bar*roe_rho[i];
	r_eigen_f[1][3]= alpha*(roe_u[i]+kx_bar*roe_a[i]);
	r_eigen_f[1][4]= alpha*(roe_u[i]-kx_bar*roe_a[i]);

	r_eigen_f[2][0]= kx_bar*roe_v[i]+kz_bar*roe_rho[i];
	r_eigen_f[2][1]= ky_bar*roe_v[i];
	r_eigen_f[2][2]= kz_bar*roe_v[i]-kx_bar*roe_rho[i];
	r_eigen_f[2][3]= alpha*(roe_v[i]+ky_bar*roe_a[i]);
	r_eigen_f[2][4]= alpha*(roe_v[i]-ky_bar*roe_a[i]);
	
	r_eigen_f[3][0]= kx_bar*roe_w[i]-ky_bar*roe_rho[i];
	r_eigen_f[3][1]= ky_bar*roe_w[i]+kx_bar*roe_rho[i];
	r_eigen_f[3][2]= kz_bar*roe_w[i];
	r_eigen_f[3][3]= alpha*(roe_w[i]+kz_bar*roe_a[i]);
	r_eigen_f[3][4]= alpha*(roe_w[i]-kz_bar*roe_a[i]);
	
	r_eigen_f[4][0]= ((kx_bar*phi_sq)/0.4)+roe_rho[i]*(kz_bar*roe_v[i]-ky_bar*roe_w[i]);
	r_eigen_f[4][1]= ((ky_bar*phi_sq)/0.4)+roe_rho[i]*(kx_bar*roe_w[i]-kz_bar*roe_u[i]);
	r_eigen_f[4][2]= ((kz_bar*phi_sq)/0.4)+roe_rho[i]*(ky_bar*roe_u[i]-kx_bar*roe_v[i]);
	r_eigen_f[4][3]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))+theta*roe_a[i]);
	r_eigen_f[4][4]= alpha*(((phi_sq+roe_a[i]*roe_a[i])/(0.4))-theta*roe_a[i]);
	
	/*********************************************/
	
	l_eigen_f[0][0]= kx_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kz_bar*roe_v[i]-ky_bar*roe_w[i])/roe_rho[i]);
	l_eigen_f[0][1]= kx_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[0][2]= (kx_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]))+(kz_bar/roe_rho[i]);
	l_eigen_f[0][3]= (kx_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))-(ky_bar/roe_rho[i]);
	l_eigen_f[0][4]= -kx_bar*0.4/(roe_a[i]*roe_a[i]);
			
	l_eigen_f[1][0]= ky_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((kx_bar*roe_w[i]-kz_bar*roe_u[i])/roe_rho[i]);
	l_eigen_f[1][1]= (ky_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i]))-kz_bar/roe_rho[i];
	l_eigen_f[1][2]= ky_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[1][3]= (ky_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]))+(kx_bar/roe_rho[i]);
	l_eigen_f[1][4]=	-ky_bar*0.4/(roe_a[i]*roe_a[i]);
				
	l_eigen_f[2][0]= kz_bar*(1.0-(phi_sq/(roe_a[i]*roe_a[i])))-((ky_bar*roe_u[i]-kx_bar*roe_v[i])/roe_rho[i]);
	l_eigen_f[2][1]= kz_bar*0.4*roe_u[i]/(roe_a[i]*roe_a[i])+(ky_bar/roe_rho[i]);
	l_eigen_f[2][2]= kz_bar*0.4*roe_v[i]/(roe_a[i]*roe_a[i])-(kx_bar/roe_rho[i]);
	l_eigen_f[2][3]= kz_bar*0.4*roe_w[i]/(roe_a[i]*roe_a[i]);
	l_eigen_f[2][4]= -kz_bar*0.4/(roe_a[i]*roe_a[i]);
	
	l_eigen_f[3][0]= beta*(phi_sq-theta*roe_a[i]);
	l_eigen_f[3][1]= -beta*(0.4*roe_u[i]-kx_bar*roe_a[i]);
	l_eigen_f[3][2]= -beta*(0.4*roe_v[i]-ky_bar*roe_a[i]);
	l_eigen_f[3][3]= -beta*(0.4*roe_w[i]-kz_bar*roe_a[i]);
	l_eigen_f[3][4]= beta*0.4;

	l_eigen_f[4][0]= beta*(phi_sq+theta*roe_a[i]);
	l_eigen_f[4][1]= -beta*(0.4*roe_u[i]+kx_bar*roe_a[i]);
	l_eigen_f[4][2]= -beta*(0.4*roe_v[i]+ky_bar*roe_a[i]);
	l_eigen_f[4][3]= -beta*(0.4*roe_w[i]+kz_bar*roe_a[i]);
	l_eigen_f[4][4]= beta*0.4;
	
	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/			

	m=0;
	for(k=0; k<5; k++)
	{	
	//	if (hk == 5)
		{
			f_iminus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][4].ip;			
			f_iminus[node[node[i].n_n[4]].n_n[4]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[4]].n_n[4]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[4]].n_n[4]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[4]].n_n[4]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[4]].n_n[4]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[4]].n_n[4]][4].ip;			
			f_iminus[node[i].n_n[4]][k] = l_eigen_f[m][0]*Q[node[i].n_n[4]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[4]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[4]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[4]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[4]][4].ip;			
			f_iminus[i][k] = l_eigen_f[m][0]*Q[i][0].ip+l_eigen_f[m][1]*Q[i][1].ip+l_eigen_f[m][2]*Q[i][2].ip+l_eigen_f[m][3]*Q[i][3].ip+l_eigen_f[m][4]*Q[i][4].ip;
			f_iminus[node[i].n_n[5]][k] = l_eigen_f[m][0]*Q[node[i].n_n[5]][0].ip+l_eigen_f[m][1]*Q[node[i].n_n[5]][1].ip+l_eigen_f[m][2]*Q[node[i].n_n[5]][2].ip+l_eigen_f[m][3]*Q[node[i].n_n[5]][3].ip+l_eigen_f[m][4]*Q[node[i].n_n[5]][4].ip;
			f_iminus[node[node[i].n_n[5]].n_n[5]][k] = l_eigen_f[m][0]*Q[node[node[i].n_n[5]].n_n[5]][0].ip+l_eigen_f[m][1]*Q[node[node[i].n_n[5]].n_n[5]][1].ip+l_eigen_f[m][2]*Q[node[node[i].n_n[5]].n_n[5]][2].ip+l_eigen_f[m][3]*Q[node[node[i].n_n[5]].n_n[5]][3].ip+l_eigen_f[m][4]*Q[node[node[i].n_n[5]].n_n[5]][4].ip;			
			f_iminus[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] = l_eigen_f[m][0]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][0].ip+l_eigen_f[m][1]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][1].ip+l_eigen_f[m][2]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][2].ip+l_eigen_f[m][3]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][3].ip+l_eigen_f[m][4]*Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][4].ip;						
		}
		
		m++; 
	}
}
	
	

	
	
	
	