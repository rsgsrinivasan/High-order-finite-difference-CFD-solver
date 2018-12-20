#include <stdio.h>
#include<stdlib.h> 
#include<math.h>
#include"function.h"
#include"mpi.h"

void roe_average(double **u, double **v, double **w, double **rho, double **p, double **t, double **mu, double **e, int g_node, ELEM *CD, JACOB *jacobian, TRANSFORMED *metric, \
int j, double **a, MNODE *node, double *roe_rho_ip, double *roe_u_ip, double *roe_v_ip, double *roe_w_ip, double *roe_h_ip, double *roe_rho_im, double *roe_u_im, double *roe_v_im, \
double *roe_w_im, double *roe_h_im, double *roe_rho_jp, double *roe_u_jp, double *roe_v_jp, double *roe_w_jp, double *roe_h_jp, double *roe_rho_jm, \
double *roe_u_jm, double *roe_v_jm, double *roe_w_jm, double *roe_h_jm, double *roe_rho_kp, double *roe_u_kp, double *roe_v_kp, double *roe_w_kp, double *roe_h_kp, \
double *roe_rho_km, double *roe_u_km, double *roe_v_km, double *roe_w_km, double *roe_h_km, int all_bou_node, int *all_boundary_nodes, RESIDUAL **div, int sd_node, int i)
{
	//int i;
	//for (i=1; i<g_node; i++)
	//{		
		
		double sqr_b1, sqr_i, sqr_b3, sqr_c0, sqr_c2, sqr_d4, sqr_d5;
		sqr_i  = sqrt(rho[j][i]);
		sqr_b1 = sqrt(rho[j][node[i].n_n[1]]);
		sqr_b3 = sqrt(rho[j][node[i].n_n[3]]);
		sqr_c0 = sqrt(rho[j][node[i].n_n[0]]);
		sqr_c2 = sqrt(rho[j][node[i].n_n[2]]);
		sqr_d4 = sqrt(rho[j][node[i].n_n[4]]);
		sqr_d5 = sqrt(rho[j][node[i].n_n[5]]);
		
		
		roe_rho_ip[i] = sqr_b1*sqr_i;
		roe_u_ip[i] = (sqr_b1*u[j][node[i].n_n[1]]+sqr_i*u[j][i])/(sqr_b1+sqr_i);
		roe_v_ip[i] = (sqr_b1*v[j][node[i].n_n[1]]+sqr_i*v[j][i])/(sqr_b1+sqr_i);
		roe_w_ip[i] = (sqr_b1*w[j][node[i].n_n[1]]+sqr_i*w[j][i])/(sqr_b1+sqr_i);
		roe_h_ip[i] = (sqr_b1*(((p[j][node[i].n_n[1]]/0.4)+0.5*rho[j][node[i].n_n[1]]*(u[j][node[i].n_n[1]]*u[j][node[i].n_n[1]]+v[j][node[i].n_n[1]]*v[j][node[i].n_n[1]]+w[j][node[i].n_n[1]]*w[j][node[i].n_n[1]])+p[j][node[i].n_n[1]])/rho[j][node[i].n_n[1]])+\
		               sqr_i*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqr_b1+sqr_i);
					
		roe_rho_im[i] = sqr_b3*sqr_i;
		roe_u_im[i] = (sqr_b3*u[j][node[i].n_n[3]]+sqr_i*u[j][i])/(sqr_b3+sqr_i);
		roe_v_im[i] = (sqr_b3*v[j][node[i].n_n[3]]+sqr_i*v[j][i])/(sqr_b3+sqr_i);
		roe_w_im[i] = (sqr_b3*w[j][node[i].n_n[3]]+sqr_i*w[j][i])/(sqr_b3+sqr_i);
		roe_h_im[i] = (sqr_b3*(((p[j][node[i].n_n[3]]/0.4)+0.5*rho[j][node[i].n_n[3]]*(u[j][node[i].n_n[3]]*u[j][node[i].n_n[3]]+v[j][node[i].n_n[3]]*v[j][node[i].n_n[3]]+w[j][node[i].n_n[3]]*w[j][node[i].n_n[3]])+p[j][node[i].n_n[3]])/rho[j][node[i].n_n[3]])+\
		               sqr_i*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqr_b3+sqr_i);
				
		roe_rho_jp[i] = sqr_c0*sqr_i;
		roe_u_jp[i] = (sqr_c0*u[j][node[i].n_n[0]]+sqr_i*u[j][i])/(sqr_c0+sqr_i);
		roe_v_jp[i] = (sqr_c0*v[j][node[i].n_n[0]]+sqr_i*v[j][i])/(sqr_c0+sqr_i);
		roe_w_jp[i] = (sqr_c0*w[j][node[i].n_n[0]]+sqr_i*w[j][i])/(sqr_c0+sqr_i);
		roe_h_jp[i] = (sqr_c0*(((p[j][node[i].n_n[0]]/0.4)+0.5*rho[j][node[i].n_n[0]]*(u[j][node[i].n_n[0]]*u[j][node[i].n_n[0]]+v[j][node[i].n_n[0]]*v[j][node[i].n_n[0]]+w[j][node[i].n_n[0]]*w[j][node[i].n_n[0]])+p[j][node[i].n_n[0]])/rho[j][node[i].n_n[0]])+\
		               sqr_i*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqr_c0+sqr_i);
					
		roe_rho_jm[i] = sqr_c2*sqr_i;
		roe_u_jm[i] = (sqr_c2*u[j][node[i].n_n[2]]+sqr_i*u[j][i])/(sqr_c2+sqr_i);
		roe_v_jm[i] = (sqr_c2*v[j][node[i].n_n[2]]+sqr_i*v[j][i])/(sqr_c2+sqr_i);
		roe_w_jm[i] = (sqr_c2*w[j][node[i].n_n[2]]+sqr_i*w[j][i])/(sqr_c2+sqr_i);
		roe_h_jm[i] = (sqr_c2*(((p[j][node[i].n_n[2]]/0.4)+0.5*rho[j][node[i].n_n[2]]*(u[j][node[i].n_n[2]]*u[j][node[i].n_n[2]]+v[j][node[i].n_n[2]]*v[j][node[i].n_n[2]]+w[j][node[i].n_n[2]]*w[j][node[i].n_n[2]])+p[j][node[i].n_n[2]])/rho[j][node[i].n_n[2]])+\
		               sqr_i*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqr_c2+sqr_i);

		roe_rho_kp[i] = sqr_d4*sqr_i;
		roe_u_kp[i] = (sqr_d4*u[j][node[i].n_n[4]]+sqr_i*u[j][i])/(sqr_d4+sqr_i);
		roe_v_kp[i] = (sqr_d4*v[j][node[i].n_n[4]]+sqr_i*v[j][i])/(sqr_d4+sqr_i);
		roe_w_kp[i] = (sqr_d4*w[j][node[i].n_n[4]]+sqr_i*w[j][i])/(sqr_d4+sqr_i);
		roe_h_kp[i] = (sqr_d4*(((p[j][node[i].n_n[4]]/0.4)+0.5*rho[j][node[i].n_n[4]]*(u[j][node[i].n_n[4]]*u[j][node[i].n_n[4]]+v[j][node[i].n_n[4]]*v[j][node[i].n_n[4]]+w[j][node[i].n_n[4]]*w[j][node[i].n_n[4]])+p[j][node[i].n_n[4]])/rho[j][node[i].n_n[4]])+\
		               sqr_i*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqr_d4+sqr_i);
					
		roe_rho_km[i] = sqr_d5*sqr_i;
		roe_u_km[i] = (sqr_d5*u[j][node[i].n_n[5]]+sqr_i*u[j][i])/(sqr_d5+sqr_i);
		roe_v_km[i] = (sqr_d5*v[j][node[i].n_n[5]]+sqr_i*v[j][i])/(sqr_d5+sqr_i);
		roe_w_km[i] = (sqr_d5*w[j][node[i].n_n[5]]+sqr_i*w[j][i])/(sqr_d5+sqr_i);
		roe_h_km[i] = (sqr_d5*(((p[j][node[i].n_n[5]]/0.4)+0.5*rho[j][node[i].n_n[5]]*(u[j][node[i].n_n[5]]*u[j][node[i].n_n[5]]+v[j][node[i].n_n[5]]*v[j][node[i].n_n[5]]+w[j][node[i].n_n[5]]*w[j][node[i].n_n[5]])+p[j][node[i].n_n[5]])/rho[j][node[i].n_n[5]])+\
		               sqr_i*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqr_d5+sqr_i);
		
		
		
/*		
		
		
		roe_rho_ip[i] = sqrt(rho[j][node[i].n_n[1]]*rho[j][i]);
		roe_u_ip[i] = (sqrt(rho[j][node[i].n_n[1]])*u[j][node[i].n_n[1]]+sqrt(rho[j][i])*u[j][i])/(sqrt(rho[j][node[i].n_n[1]])+sqrt(rho[j][i]));
		roe_v_ip[i] = (sqrt(rho[j][node[i].n_n[1]])*v[j][node[i].n_n[1]]+sqrt(rho[j][i])*v[j][i])/(sqrt(rho[j][node[i].n_n[1]])+sqrt(rho[j][i]));
		roe_w_ip[i] = (sqrt(rho[j][node[i].n_n[1]])*w[j][node[i].n_n[1]]+sqrt(rho[j][i])*w[j][i])/(sqrt(rho[j][node[i].n_n[1]])+sqrt(rho[j][i]));
		roe_h_ip[i] = (sqrt(rho[j][node[i].n_n[1]])*(((p[j][node[i].n_n[1]]/0.4)+0.5*rho[j][node[i].n_n[1]]*(u[j][node[i].n_n[1]]*u[j][node[i].n_n[1]]+v[j][node[i].n_n[1]]*v[j][node[i].n_n[1]]+w[j][node[i].n_n[1]]*w[j][node[i].n_n[1]])+p[j][node[i].n_n[1]])/rho[j][node[i].n_n[1]])+\
		               sqrt(rho[j][i])*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqrt(rho[j][node[i].n_n[1]])+sqrt(rho[j][i]));
					
		roe_rho_im[i] = sqrt(rho[j][node[i].n_n[3]]*rho[j][i]);
		roe_u_im[i] = (sqrt(rho[j][node[i].n_n[3]])*u[j][node[i].n_n[3]]+sqrt(rho[j][i])*u[j][i])/(sqrt(rho[j][node[i].n_n[3]])+sqrt(rho[j][i]));
		roe_v_im[i] = (sqrt(rho[j][node[i].n_n[3]])*v[j][node[i].n_n[3]]+sqrt(rho[j][i])*v[j][i])/(sqrt(rho[j][node[i].n_n[3]])+sqrt(rho[j][i]));
		roe_w_im[i] = (sqrt(rho[j][node[i].n_n[3]])*w[j][node[i].n_n[3]]+sqrt(rho[j][i])*w[j][i])/(sqrt(rho[j][node[i].n_n[3]])+sqrt(rho[j][i]));
		roe_h_im[i] = (sqrt(rho[j][node[i].n_n[3]])*(((p[j][node[i].n_n[3]]/0.4)+0.5*rho[j][node[i].n_n[3]]*(u[j][node[i].n_n[3]]*u[j][node[i].n_n[3]]+v[j][node[i].n_n[3]]*v[j][node[i].n_n[3]]+w[j][node[i].n_n[3]]*w[j][node[i].n_n[3]])+p[j][node[i].n_n[3]])/rho[j][node[i].n_n[3]])+\
		               sqrt(rho[j][i])*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqrt(rho[j][node[i].n_n[3]])+sqrt(rho[j][i]));
				
		roe_rho_jp[i] = sqrt(rho[j][node[i].n_n[0]]*rho[j][i]);
		roe_u_jp[i] = (sqrt(rho[j][node[i].n_n[0]])*u[j][node[i].n_n[0]]+sqrt(rho[j][i])*u[j][i])/(sqrt(rho[j][node[i].n_n[0]])+sqrt(rho[j][i]));
		roe_v_jp[i] = (sqrt(rho[j][node[i].n_n[0]])*v[j][node[i].n_n[0]]+sqrt(rho[j][i])*v[j][i])/(sqrt(rho[j][node[i].n_n[0]])+sqrt(rho[j][i]));
		roe_w_jp[i] = (sqrt(rho[j][node[i].n_n[0]])*w[j][node[i].n_n[0]]+sqrt(rho[j][i])*w[j][i])/(sqrt(rho[j][node[i].n_n[0]])+sqrt(rho[j][i]));
		roe_h_jp[i] = (sqrt(rho[j][node[i].n_n[0]])*(((p[j][node[i].n_n[0]]/0.4)+0.5*rho[j][node[i].n_n[0]]*(u[j][node[i].n_n[0]]*u[j][node[i].n_n[0]]+v[j][node[i].n_n[0]]*v[j][node[i].n_n[0]]+w[j][node[i].n_n[0]]*w[j][node[i].n_n[0]])+p[j][node[i].n_n[0]])/rho[j][node[i].n_n[0]])+\
		               sqrt(rho[j][i])*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqrt(rho[j][node[i].n_n[0]])+sqrt(rho[j][i]));
					
		roe_rho_jm[i] = sqrt(rho[j][node[i].n_n[2]]*rho[j][i]);
		roe_u_jm[i] = (sqrt(rho[j][node[i].n_n[2]])*u[j][node[i].n_n[2]]+sqrt(rho[j][i])*u[j][i])/(sqrt(rho[j][node[i].n_n[2]])+sqrt(rho[j][i]));
		roe_v_jm[i] = (sqrt(rho[j][node[i].n_n[2]])*v[j][node[i].n_n[2]]+sqrt(rho[j][i])*v[j][i])/(sqrt(rho[j][node[i].n_n[2]])+sqrt(rho[j][i]));
		roe_w_jm[i] = (sqrt(rho[j][node[i].n_n[2]])*w[j][node[i].n_n[2]]+sqrt(rho[j][i])*w[j][i])/(sqrt(rho[j][node[i].n_n[2]])+sqrt(rho[j][i]));
		roe_h_jm[i] = (sqrt(rho[j][node[i].n_n[2]])*(((p[j][node[i].n_n[2]]/0.4)+0.5*rho[j][node[i].n_n[2]]*(u[j][node[i].n_n[2]]*u[j][node[i].n_n[2]]+v[j][node[i].n_n[2]]*v[j][node[i].n_n[2]]+w[j][node[i].n_n[2]]*w[j][node[i].n_n[2]])+p[j][node[i].n_n[2]])/rho[j][node[i].n_n[2]])+\
		               sqrt(rho[j][i])*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqrt(rho[j][node[i].n_n[2]])+sqrt(rho[j][i]));

		roe_rho_kp[i] = sqrt(rho[j][node[i].n_n[4]]*rho[j][i]);
		roe_u_kp[i] = (sqrt(rho[j][node[i].n_n[4]])*u[j][node[i].n_n[4]]+sqrt(rho[j][i])*u[j][i])/(sqrt(rho[j][node[i].n_n[4]])+sqrt(rho[j][i]));
		roe_v_kp[i] = (sqrt(rho[j][node[i].n_n[4]])*v[j][node[i].n_n[4]]+sqrt(rho[j][i])*v[j][i])/(sqrt(rho[j][node[i].n_n[4]])+sqrt(rho[j][i]));
		roe_w_kp[i] = (sqrt(rho[j][node[i].n_n[4]])*w[j][node[i].n_n[4]]+sqrt(rho[j][i])*w[j][i])/(sqrt(rho[j][node[i].n_n[4]])+sqrt(rho[j][i]));
		roe_h_kp[i] = (sqrt(rho[j][node[i].n_n[4]])*(((p[j][node[i].n_n[4]]/0.4)+0.5*rho[j][node[i].n_n[4]]*(u[j][node[i].n_n[4]]*u[j][node[i].n_n[4]]+v[j][node[i].n_n[4]]*v[j][node[i].n_n[4]]+w[j][node[i].n_n[4]]*w[j][node[i].n_n[4]])+p[j][node[i].n_n[4]])/rho[j][node[i].n_n[4]])+\
		               sqrt(rho[j][i])*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqrt(rho[j][node[i].n_n[4]])+sqrt(rho[j][i]));
					
		roe_rho_km[i] = sqrt(rho[j][node[i].n_n[5]]*rho[j][i]);
		roe_u_km[i] = (sqrt(rho[j][node[i].n_n[5]])*u[j][node[i].n_n[5]]+sqrt(rho[j][i])*u[j][i])/(sqrt(rho[j][node[i].n_n[5]])+sqrt(rho[j][i]));
		roe_v_km[i] = (sqrt(rho[j][node[i].n_n[5]])*v[j][node[i].n_n[5]]+sqrt(rho[j][i])*v[j][i])/(sqrt(rho[j][node[i].n_n[5]])+sqrt(rho[j][i]));
		roe_w_km[i] = (sqrt(rho[j][node[i].n_n[5]])*w[j][node[i].n_n[5]]+sqrt(rho[j][i])*w[j][i])/(sqrt(rho[j][node[i].n_n[5]])+sqrt(rho[j][i]));
		roe_h_km[i] = (sqrt(rho[j][node[i].n_n[5]])*(((p[j][node[i].n_n[5]]/0.4)+0.5*rho[j][node[i].n_n[5]]*(u[j][node[i].n_n[5]]*u[j][node[i].n_n[5]]+v[j][node[i].n_n[5]]*v[j][node[i].n_n[5]]+w[j][node[i].n_n[5]]*w[j][node[i].n_n[5]])+p[j][node[i].n_n[5]])/rho[j][node[i].n_n[5]])+\
		               sqrt(rho[j][i])*(((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i])+p[j][i])/rho[j][i]))/(sqrt(rho[j][node[i].n_n[5]])+sqrt(rho[j][i]));
	*/	
		if(i <= sd_node)
		{               
			div[i][0].u = u[0][i];
			div[i][0].v = v[0][i];
			div[i][0].e = e[0][i]; 
		}
	//}
}
