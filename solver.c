#include <stdio.h>
#include<stdlib.h> 
#include<math.h>
#include<time.h>
#include"function.h"
#include"mpi.h"


void solver(ELEM *CD, JACOB *jacobian, TRANSFORMED *metric, double ***U, double **F1, double **E1, double **G1, double **Fv1, double **Ev1, double **Gv1, double **E, double **F, double **G, double **Ev, double **Fv, double **Gv, \
double *roe_rho_ip, double *roe_u_ip, double *roe_v_ip, double *roe_w_ip, double *roe_h_ip, double *roe_rho_im, double *roe_u_im, double *roe_v_im, double *roe_w_im, double *roe_h_im, double *roe_rho_jp, double *roe_u_jp, \
double *roe_v_jp, double *roe_w_jp, double *roe_h_jp, double *roe_rho_jm, double *roe_u_jm, double *roe_v_jm, double *roe_w_jm, double *roe_h_jm, double *roe_rho_kp, double *roe_u_kp, double *roe_v_kp, double *roe_w_kp, \
double *roe_h_kp, double *roe_rho_km, double *roe_u_km, double *roe_v_km, double *roe_w_km, double *roe_h_km, double **r_eigen_Qip, double **r_eigen_Qjp, double **r_eigen_Qkp, double **r_eigen_Qim, double **r_eigen_Qjm, \
double **r_eigen_Qkm, double **l_eigen_Qip, double **l_eigen_Qjp, double **l_eigen_Qkp, double **Qi_iplus, double **Qi_iminus, double **Qj_iplus, double **Qj_iminus, double **Qk_iplus,\
double **Qk_iminus, double **Qi_half_p, double **Qi_half_np, double **Qi_half_m, double **Qi_halfn_m, double **Qj_half_p, double **Qj_half_np, double **Qj_half_m, double **Qj_halfn_m, double **Qk_half_p, \
double **Qk_half_np, double **Qk_half_m, double **Qk_halfn_m, double **IS_Qim, double **IS_Qinm, double **IS_Qjm, double **IS_Qjnm,  double **IS_Qkm, double **IS_Qknm, double **IS_Qip, double **IS_Qinp, \
double **IS_Qjp, double **IS_Qjnp, double **IS_Qkp, double **IS_Qknp, double **w_Qip, double **w_Qinp, double **w_Qim, double **w_Qinm, double **w_Qjp, double **w_Qjnp, double **w_Qjm, double **w_Qjnm,  \
double **w_Qkp, double **w_Qknp, double **w_Qkm, double **w_Qknm, double **W_Qip, double **W_Qinp, double **W_Qim, double **W_Qinm, double **W_Qjp, double **W_Qjnp, double **W_Qjm, double **W_Qjnm, \
double **W_Qkp, double **W_Qknp, double **W_Qkm, double **W_Qknm, double *Qi_iplus_half, double *Qi_iminus_half, double *Qj_iplus_half, double *Qj_iminus_half,  double *Qk_iplus_half, double *Qk_iminus_half, \
double *dF, double *dE, double *dG, double *dFv, double *dEv, double *dGv, double ***L, int NUMNP, int NELEM, double **rho, double **u, double **v, double **w, double **e, double **p, double **t, \
double **a, double **mu, int g_node, int g_elem, int cor, int all_bou_node, int *inlet, int *outlet, int *wall, int *boundary, int iterations, double *tauzz,  double *tauee,  double *tauxx, double *tauze, double *tauzx, double *tauex, double *qz, double *qe, \
double *qx, double *det, int *all_boundary_nodes, double Reyl, double **Qi_iplus_n, double **Qi_iminus_n, double **Qj_iplus_n, double **Qj_iminus_n, double **Qk_iplus_n, double **Qk_iminus_n, double *final_U, \
int inl, int out, int wal, int bou, int *inlet_node, int *outlet_node, int *wall_node, int *boundary_node, int inl_elem, int out_elem, int bou_elem, int wal_elem, int nodes, int wal_node, MNODE *node, double *roe_R, \
double *roe_a, int out_node, int inl_node, int bou_node, double *v_dash1, double *del_cfl, RESIDUAL **div, double **l_eigen_Qim, double **l_eigen_Qjm, double **l_eigen_Qkm, int sd_node, int **loc_dat, int **recv_b, int *c, \
int *proc_node, int neigh_pro, int *recv_c, int itera, int no_of_tip_send, int no_of_tip_recv, DOM_TIP *tip, DOM_TIP *tip_recv, int sd_wal_node,int sd_bou_node, int sd_out_node, int sd_inl_node, int *sd_wall_node, \
int *sd_boundary_node, int *sd_outlet_node, int *sd_inlet_node, singul *singular, Qip **Q, DETERM *deter, double *TAU_SGS_XY, double *TAU_SGS_YZ, double *TAU_SGS_XZ,\
double *H_SGS_X, double *H_SGS_Y, double *H_SGS_Z, double *D_SGS_X, double *D_SGS_Y, double *D_SGS_Z, double *TAU_SGS_XX, double *TAU_SGS_YY, double *TAU_SGS_ZZ, double *DUCROS)
{
	int i,j,k,m, iter, weno_exec, initial, gar,lm, value, u_loc, v_loc, e_loc; 
	time_t start,end;
	double epsi, temp, rc;
	//double del_fp, del_fm;
	int position, memsize, myrank, size, memsize1;
	char *buffer, *buffer2;	
	
	double Qi_iplus_half_pos[6], Qi_iminus_half_pos[6], Qi_iplus_half_neg[6], Qi_iminus_half_neg[6], Qi_iplus_half_f[6], Qi_iminus_half_f[6];
	double Qj_iplus_half_pos[6], Qj_iminus_half_pos[6], Qj_iplus_half_neg[6], Qj_iminus_half_neg[6], Qj_iplus_half_E[6], Qj_iminus_half_E[6];
	double Qk_iplus_half_pos[6], Qk_iminus_half_pos[6], Qk_iplus_half_neg[6], Qk_iminus_half_neg[6], Qk_iplus_half_G[6], Qk_iminus_half_G[6];
	double Qi_iplus_half_pos_char[6], Qi_iminus_half_pos_char[6], Qi_iplus_half_neg_char[6], Qi_iminus_half_neg_char[6];
	double Qj_iplus_half_pos_char[6], Qj_iminus_half_pos_char[6], Qj_iplus_half_neg_char[6], Qj_iminus_half_neg_char[6];
	double Qk_iplus_half_pos_char[6], Qk_iminus_half_pos_char[6], Qk_iplus_half_neg_char[6], Qk_iminus_half_neg_char[6];
	double h_F_ip[6], h_F_im[6], h_E_ip[6], h_E_im[6], h_G_ip[6], h_G_im[6], d2F_d2z_ip[6], d2F_d2z_im[6], d2E_d2e_ip[6], d2E_d2e_im[6], d2G_d2x_ip[6], d2G_d2x_im[6];  
	double F_ip_h[6], F_im_h[6], E_ip_h[6], E_im_h[6], G_ip_h[6], G_im_h[6], d4F_d4z_ip[6], d4F_d4z_im[6], d4E_d4e_ip[6], d4E_d4e_im[6], d4G_d4x_ip[6], d4G_d4x_im[6]; 
	double F_ip_pos[6], F_ip_neg[6], F_im_pos[6], F_im_neg[6], E_ip_pos[6], E_ip_neg[6], E_im_pos[6], E_im_neg[6], G_ip_pos[6], G_ip_neg[6], G_im_pos[6], G_im_neg[6];
	double F_jp_pos[6], F_jp_neg[6], F_jm_pos[6], F_jm_neg[6], E_jp_pos[6], E_jp_neg[6], E_jm_pos[6], E_jm_neg[6], G_jp_pos[6], G_jp_neg[6], G_jm_pos[6], G_jm_neg[6];
	double F_kp_pos[6], F_kp_neg[6], F_km_pos[6], F_km_neg[6], E_kp_pos[6], E_kp_neg[6], E_km_pos[6], E_km_neg[6], G_kp_pos[6], G_kp_neg[6], G_km_pos[6], G_km_neg[6];
	double F_ip_pos_comp[6], F_ip_neg_comp[6], F_im_pos_comp[6], F_im_neg_comp[6], E_ip_pos_comp[6], E_ip_neg_comp[6], E_im_pos_comp[6], E_im_neg_comp[6], G_ip_pos_comp[6], G_ip_neg_comp[6], G_im_pos_comp[6], G_im_neg_comp[6];
	double F_jp_pos_comp[6], F_jp_neg_comp[6], F_jm_pos_comp[6], F_jm_neg_comp[6], E_jp_pos_comp[6], E_jp_neg_comp[6], E_jm_pos_comp[6], E_jm_neg_comp[6], G_jp_pos_comp[6], G_jp_neg_comp[6], G_jm_pos_comp[6], G_jm_neg_comp[6];
	double F_kp_pos_comp[6], F_kp_neg_comp[6], F_km_pos_comp[6], F_km_neg_comp[6], E_kp_pos_comp[6], E_kp_neg_comp[6], E_km_pos_comp[6], E_km_neg_comp[6], G_kp_pos_comp[6], G_kp_neg_comp[6], G_km_pos_comp[6], G_km_neg_comp[6];
	double F_ip[6], F_im[6], E_ip[6], E_im[6], G_ip[6], G_im[6], F_jp[6], F_jm[6], E_jp[6], E_jm[6], G_jp[6], G_jm[6], F_kp[6], F_km[6], E_kp[6], E_km[6], G_kp[6], G_km[6];
	double *alpha_u_ip, *alpha_v_ip,  *alpha_w_ip, *alpha_u_im, *alpha_v_im, *alpha_w_im;
	double dF_W[6], dF_C[6], dE_W[6], dE_C[6], dG_W[6], dG_C[6];
	//double f_dash_Qip[8], f_dash_Qim[8], E_dash_Qjp[8], E_dash_Qjm[8], G_dash_Qkp[8], G_dash_Qkm[8], max_u, max_v, max_w;
	double *alpha_u_jp, *alpha_v_jp, *alpha_w_jp, *alpha_u_jm, *alpha_v_jm, *alpha_w_jm, *alpha_u_kp, *alpha_v_kp, *alpha_w_kp, *alpha_u_km, *alpha_v_km, *alpha_w_km, *alpha_u, *alpha_v, *alpha_w;
	//double f_dash_Qinp[8], f_dash_Qinm[8], E_dash_Qjnp[8], E_dash_Qjnm[8], Q_dash_Qknp[8], G_dash_Qknm[8];
	double eigen_Qip[8], eigen_Qim[8], eigen_Qjp[8], eigen_Qjm[8], eigen_Qkp[8], eigen_Qkm[8], eigen_Qinp[8], eigen_Qinm[8], eigen_Qjnp[8], eigen_Qjnm[8], eigen_Qknp[8], eigen_Qknm[8];
	//double eigen_QEip[8], eigen_QEim[8], eigen_Qfjp[8], eigen_Qfjm[8], eigen_QGkp[8], eigen_QGkm[8], eigen_QEinp[8], eigen_QEinm[8], eigen_Qfjnp[8], eigen_Qfjnm[8], eigen_QGknp[8], eigen_QGknm[8];
	//double rj_F, rj_E, rjp_F, rjp_E, rj_p_F, rj_p_E, epsi_F, epsi_E, del_F_iphalf, del_F_imhalf, del_F_ip3half, del_E_iphalf, del_E_imhalf, del_E_ip3half, zeeta, eeta;
	//double rj_G, rjp_G, rj_p_G, epsi_G, del_G_iphalf, del_G_imhalf, del_G_ip3half;
	
	double zeta_xip, zeta_xim, zeta_yip, zeta_yim, zeta_zip, zeta_zim, eta_xjp, eta_xjm, eta_yjp, eta_yjm, eta_zjp, eta_zjm, xi_xkp, xi_xkm, xi_ykp, xi_ykm, xi_zkp, xi_zkm; 
	double e_inf, max_div_u, max_div_v, max_div_w, max_div_e;
//	double *flux_det_ip, *flux_det_im, *flux_det_jp, *flux_det_jm, *flux_det_kp, *flux_det_km;
	e_inf = e[0][10650];
	epsi = 0.003;
	FILE *fnode;
	char line[1000];
	char filename[100];
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Comm_size(MPI_COMM_WORLD,&size);
  	MPI_Status status;
	rc = 0.3;
//	e_inf = e[0][10650];
	epsi = 0.003;
	
	max_div_u =-1e25;
	max_div_v =-1e25;
	max_div_w =-1e25;
	max_div_e =-1e25;
	
	
/*	alpha_u = (double **)malloc(g_node*sizeof(double*));
	alpha_v = (double **)malloc(g_node*sizeof(double*));
	alpha_w = (double **)malloc(g_node*sizeof(double*));
	alpha_u_ip = (double **)malloc(g_node*sizeof(double*));
	alpha_u_im = (double **)malloc(g_node*sizeof(double*));
	alpha_v_ip = (double **)malloc(g_node*sizeof(double*));
	alpha_v_im = (double **)malloc(g_node*sizeof(double*));
	alpha_w_ip = (double **)malloc(g_node*sizeof(double*));
	alpha_w_im = (double **)malloc(g_node*sizeof(double*));
	alpha_u_jp = (double **)malloc(g_node*sizeof(double*));
	alpha_u_jm = (double **)malloc(g_node*sizeof(double*));
	alpha_v_jp = (double **)malloc(g_node*sizeof(double*));
	alpha_v_jm = (double **)malloc(g_node*sizeof(double*));
	alpha_w_jp = (double **)malloc(g_node*sizeof(double*));
	alpha_w_jm = (double **)malloc(g_node*sizeof(double*));
	alpha_u_kp = (double **)malloc(g_node*sizeof(double*));
	alpha_u_km = (double **)malloc(g_node*sizeof(double*));
	alpha_v_kp = (double **)malloc(g_node*sizeof(double*));
	alpha_v_km = (double **)malloc(g_node*sizeof(double*));
	alpha_w_kp = (double **)malloc(g_node*sizeof(double*));
	alpha_w_km = (double **)malloc(g_node*sizeof(double*));
	flux_det_ip = (double *)malloc(g_node*sizeof(double));
	flux_det_im = (double *)malloc(g_node*sizeof(double));
	flux_det_jp = (double *)malloc(g_node*sizeof(double));
	flux_det_jm = (double *)malloc(g_node*sizeof(double));
	flux_det_kp = (double *)malloc(g_node*sizeof(double));
	flux_det_km = (double *)malloc(g_node*sizeof(double));
	*/
	
		
	alpha_u = (double *)malloc(5*sizeof(double));
	alpha_v = (double *)malloc(5*sizeof(double));
	alpha_w = (double *)malloc(5*sizeof(double));
	alpha_u_ip = (double *)malloc(5*sizeof(double));
	alpha_u_im = (double *)malloc(5*sizeof(double));
	alpha_v_ip = (double *)malloc(5*sizeof(double));
	alpha_v_im = (double *)malloc(5*sizeof(double));
	alpha_w_ip = (double *)malloc(5*sizeof(double));
	alpha_w_im = (double *)malloc(5*sizeof(double));
	alpha_u_jp = (double *)malloc(5*sizeof(double));
	alpha_u_jm = (double *)malloc(5*sizeof(double));
	alpha_v_jp = (double *)malloc(5*sizeof(double));
	alpha_v_jm= (double *)malloc(5*sizeof(double));
	alpha_w_jp = (double *)malloc(5*sizeof(double));
	alpha_w_jm = (double *)malloc(5*sizeof(double));
	alpha_u_kp = (double *)malloc(5*sizeof(double));
	alpha_u_km = (double *)malloc(5*sizeof(double));
	alpha_v_kp = (double *)malloc(5*sizeof(double));
	alpha_v_km = (double *)malloc(5*sizeof(double));
	alpha_w_kp = (double *)malloc(5*sizeof(double));
	alpha_w_km = (double *)malloc(5*sizeof(double));

	
	lm = 0;
	for (iter=itera; iter<iterations; iter++)
	{	
		start = time(NULL);	
		max_div_u =-1e25;
		max_div_v =-1e25;
		max_div_e =-1e25;
		gar =0;
		
		weno_exec = 1;
		
		for(j=0; j<3; j++)
		{
			lm++;
			gar++;
			for (i=1; i<g_node; i++)
			{	
				roe_average(u, v, w, rho, p, t, mu, e, g_node, CD, jacobian, metric, j, a, node, roe_rho_ip, roe_u_ip, roe_v_ip, roe_w_ip, roe_h_ip, roe_rho_im, roe_u_im, roe_v_im, \
				roe_w_im, roe_h_im, roe_rho_jp, roe_u_jp, roe_v_jp, roe_w_jp, roe_h_jp, roe_rho_jm, roe_u_jm, roe_v_jm, roe_w_jm, roe_h_jm, \
				roe_rho_kp, roe_u_kp, roe_v_kp, roe_w_kp, roe_h_kp, roe_rho_km, roe_u_km, roe_v_km, roe_w_km, roe_h_km, all_bou_node, all_boundary_nodes, div, sd_node, i);
			
				U[i][j][0] = (det[i])*rho[j][i];
				U[i][j][1] = (det[i])*rho[j][i]*u[j][i];
				U[i][j][2] = (det[i])*rho[j][i]*v[j][i];
				U[i][j][3] = (det[i])*rho[j][i]*w[j][i];
				U[i][j][4] = (det[i])*((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i]));
				
				Q[i][0].ip = rho[j][i];
				Q[i][1].ip = rho[j][i]*u[j][i];
				Q[i][2].ip = rho[j][i]*v[j][i];
				Q[i][3].ip = rho[j][i]*w[j][i];
				Q[i][4].ip = ((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i]));
			}
		
			for (i=1; i<g_node; i++)
			{			
				viscousflux_variables(CD, u, v, w, t, j, NUMNP, node, mu, Reyl, tauzz, tauze, tauee, tauxx, tauzx, tauex, qz, qe, qx, metric, g_node, i, wall_node, U, singular, TAU_SGS_XY, TAU_SGS_YZ, TAU_SGS_XZ,\
				H_SGS_X, H_SGS_Y, H_SGS_Z, D_SGS_X, D_SGS_Y, D_SGS_Z, TAU_SGS_XX, TAU_SGS_YY, TAU_SGS_ZZ, det, DUCROS);
				
				F1[i][0] = rho[j][i]*u[j][i];
				F1[i][1] = rho[j][i]*u[j][i]*u[j][i]+p[j][i];
				F1[i][2] = rho[j][i]*v[j][i]*u[j][i];
				F1[i][3] = rho[j][i]*u[j][i]*w[j][i];
				F1[i][4] = (((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i]))+p[j][i])*u[j][i];
		
				E1[i][0] = rho[j][i]*v[j][i];
				E1[i][1] = rho[j][i]*v[j][i]*u[j][i];
				E1[i][2] = rho[j][i]*v[j][i]*v[j][i]+p[j][i];
				E1[i][3] = rho[j][i]*v[j][i]*w[j][i];
				E1[i][4] = (((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i]))+p[j][i])*v[j][i];
		
				G1[i][0] = rho[j][i]*w[j][i];
				G1[i][1] = rho[j][i]*u[j][i]*w[j][i];
				G1[i][2] = rho[j][i]*v[j][i]*w[j][i];
				G1[i][3] = rho[j][i]*w[j][i]*w[j][i]+p[j][i];
				G1[i][4] = (((p[j][i]/0.4)+0.5*rho[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]+w[j][i]*w[j][i]))+p[j][i])*w[j][i];

				Fv1[i][0] = 0.0;
				Fv1[i][1] = (-1.0)*tauzz[i]+TAU_SGS_XX[i];
				Fv1[i][2] = (-1.0)*tauze[i]+TAU_SGS_XY[i];
				Fv1[i][3] = (-1.0)*tauzx[i]+TAU_SGS_XZ[i];
				Fv1[i][4] = (((-1.0)*u[j][i]*tauzz[i]-v[j][i]*tauze[i]-w[j][i]*tauzx[i])+qz[i]+H_SGS_X[i]+D_SGS_X[i]);
		
				Ev1[i][0] = 0.0;
				Ev1[i][1] = (-1.0)*tauze[i]+TAU_SGS_XY[i];
				Ev1[i][2] = (-1.0)*tauee[i]+TAU_SGS_YY[i];
				Ev1[i][3] = (-1.0)*tauex[i]+TAU_SGS_YZ[i];
				Ev1[i][4] = (((-1.0)*u[j][i]*tauze[i]-v[j][i]*tauee[i]-w[j][i]*tauex[i])+qe[i]+H_SGS_Y[i]+D_SGS_Y[i]);
				
				Gv1[i][0] = 0.0;
				Gv1[i][1] = (-1.0)*tauzx[i]+TAU_SGS_XZ[i];
				Gv1[i][2] = (-1.0)*tauex[i]+TAU_SGS_YZ[i];
				Gv1[i][3] = (-1.0)*tauxx[i]+TAU_SGS_ZZ[i];
				Gv1[i][4] = (((-1.0)*u[j][i]*tauzx[i]-v[j][i]*tauex[i]-w[j][i]*tauxx[i])+qx[i]+H_SGS_Z[i]+D_SGS_Z[i]);
				
				for(k=0; k<5; k++)
				{
					F[i][k] =(det[i])*(F1[i][k]*metric[i].zeta_x+E1[i][k]*metric[i].zeta_y+G1[i][k]*metric[i].zeta_z);
					E[i][k] =(det[i])*(F1[i][k]*metric[i].eta_x+E1[i][k]*metric[i].eta_y+G1[i][k]*metric[i].eta_z);
					G[i][k] =(det[i])*(F1[i][k]*metric[i].xi_x+E1[i][k]*metric[i].xi_y+G1[i][k]*metric[i].xi_z);
					
					Fv[i][k] =(det[i])*(Fv1[i][k]*metric[i].zeta_x+Ev1[i][k]*metric[i].zeta_y+Gv1[i][k]*metric[i].zeta_z);
					Ev[i][k] =(det[i])*(Fv1[i][k]*metric[i].eta_x+Ev1[i][k]*metric[i].eta_y+Gv1[i][k]*metric[i].eta_z);
					Gv[i][k] =(det[i])*(Fv1[i][k]*metric[i].xi_x+Ev1[i][k]*metric[i].xi_y+Gv1[i][k]*metric[i].xi_z);
				}
				node[i].val = 0;
			}	
/*			
			sprintf(filename,"Fv_%d.txt",myrank);
			fnode=fopen(filename,"w");
			for (i=1; i<=sd_node; i++)
			{
				fprintf(fnode,"%e %e %e %e %e %e %d %d\n", tauzz[i], tauze[i], tauzx[i], tauze[i], tauee[i], tauex[i], node[i].global, node[i].loc);
			}		
			fclose(fnode);
			
			sprintf(filename,"Ev_%d.txt",myrank);
			fnode=fopen(filename,"w");
			for (i=1; i<=sd_node; i++)
			{
				fprintf(fnode,"%e %e %e %e %e %e %e %d %d\n", tauzx[i], tauex[i], Ev[i][2], tauxx[i], qz[i], qe[i], qx[i], node[i].global, node[i].loc);
			}		
			fclose(fnode);
			
			sprintf(filename,"Gv_%d.txt",myrank);
			fnode=fopen(filename,"w");
			for (i=1; i<=sd_node; i++)
			{
				fprintf(fnode,"%e %e %e %e %e %d %d\n", Gv[i][0], Gv[i][1], Gv[i][2], Gv[i][3], U[i][j][4], node[i].global, node[i].loc);
			}		
			fclose(fnode);
			
			sprintf(filename,"deriv_%d.txt",myrank);
			fnode=fopen(filename,"w");
			
	*/		
/********************************************************************************************************************************************************************************************/			
			for (i=1; i<=sd_node; i++)
			{
				
			/*	if (DUCROS[node[i].n_n[1]] >=0.65 && node[i].val == 0)
				{
					DUCROS[i] = DUCROS[node[i].n_n[1]];
					node[i].val = 1;
				}
				
				if (DUCROS[node[node[i].n_n[1]].n_n[1]] >=0.65 && node[i].val == 0)
				{
					DUCROS[i] = DUCROS[node[node[i].n_n[1]].n_n[1]];
					node[i].val = 1;
				}
				
				if (DUCROS[node[i].n_n[3]] >=0.65 && node[i].val == 0)
				{
					DUCROS[i] = DUCROS[node[i].n_n[3]];
					node[i].val = 1;
				}
				
				if (DUCROS[node[node[i].n_n[3]].n_n[3]] >=0.65 && node[i].val == 0)
				{
					DUCROS[i] = DUCROS[node[node[i].n_n[3]].n_n[3]];
					node[i].val = 1;
				}
				
			*/	
			//	if (node[i].loc == 0 && node[i].corner_ID == 0 )
				{
					value = node[i].val;
					
					weno_solver_1(CD, roe_u_ip, roe_v_ip, roe_w_ip, roe_h_ip, roe_rho_ip, r_eigen_Qip, l_eigen_Qip, Q, j, 1, g_elem, u, v, w, e, p, rho, U,  \
					g_node, NUMNP, NELEM, roe_R, roe_a, node, i, metric, deter, Qi_iplus, Qi_iminus_n);
					
					weno_solver_3(CD, roe_u_im, roe_v_im, roe_w_im, roe_h_im, roe_rho_im, r_eigen_Qim, l_eigen_Qim, Q, j, 3, g_elem, u, v, w, e, p, rho, U,  \
					g_node, NUMNP, NELEM, roe_R, roe_a, node, i, metric, deter, Qi_iplus, Qi_iminus_n);
					
					weno_solver_0(CD, roe_u_jp, roe_v_jp, roe_w_jp, roe_h_jp, roe_rho_jp, r_eigen_Qjp, l_eigen_Qjp, Q, j, 0, g_elem, u, v, w, e, p, rho, U,  \
					g_node, NUMNP, NELEM, roe_R, roe_a, node, i, metric, deter, Qj_iplus, Qj_iminus_n);
					
					weno_solver_2(CD, roe_u_jm, roe_v_jm, roe_w_jm, roe_h_jm, roe_rho_jm, r_eigen_Qjm, l_eigen_Qjm, Q, j, 2, g_elem, u, v, w, e, p, rho, U,  \
					g_node, NUMNP, NELEM, roe_R, roe_a, node, i, metric, deter, Qj_iplus, Qj_iminus_n);
					
					weno_solver_4(CD, roe_u_kp, roe_v_kp, roe_w_kp, roe_h_kp, roe_rho_kp, r_eigen_Qkp, l_eigen_Qkp, Q, j, 4, g_elem, u, v, w, e, p, rho, U,  \
					g_node, NUMNP, NELEM, roe_R, roe_a, node, i, metric, deter, Qk_iplus, Qk_iminus_n);
					
					weno_solver_5(CD, roe_u_km, roe_v_km, roe_w_km, roe_h_km, roe_rho_km, r_eigen_Qkm, l_eigen_Qkm, Q, j, 5, g_elem, u, v, w, e, p, rho, U,  \
					g_node, NUMNP, NELEM, roe_R, roe_a, node, i, metric, deter, Qk_iplus, Qk_iminus_n);
				
					for (k=0;k<5;k++)
					{		
						/** Qi(i+1/2)_plus**/
						Qi_half_p[0][k] = (15.0/8.0)*Qi_iplus[node[i].n_n[1]][k]-(5.0/4.0)*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k]+(3.0/8.0)*Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k];           
						Qi_half_p[1][k] = (3.0/8.0)*Qi_iplus[i][k]+(3.0/4.0)*Qi_iplus[node[i].n_n[1]][k]-(1.0/8.0)*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k];
						Qi_half_p[2][k] = (-1.0/8.0)*Qi_iplus[node[i].n_n[3]][k]+(3.0/4.0)*Qi_iplus[i][k]+(3.0/8.0)*Qi_iplus[node[i].n_n[1]][k];
						
						/** Qi(i+1/2)_minus**/
						Qi_half_m[0][k] = (3.0/8.0)*Qi_iplus[i][k]+(3.0/4.0)*Qi_iplus[node[i].n_n[1]][k]-(1.0/8.0)*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k];  
						Qi_half_m[1][k] = (-1.0/8.0)*Qi_iplus[node[i].n_n[3]][k]+(3.0/4.0)*Qi_iplus[i][k]+(3.0/8.0)*Qi_iplus[node[i].n_n[1]][k];
						Qi_half_m[2][k] = (3.0/8.0)*Qi_iplus[node[node[i].n_n[3]].n_n[3]][k]-(5.0/4.0)*Qi_iplus[node[i].n_n[3]][k]+(15.0/8.0)*Qi_iplus[i][k];
						
						/** Qi(i-1/2)_minus**/
						Qi_halfn_m[0][k] = (3.0/8.0)*Qi_iminus_n[node[i].n_n[3]][k]+(3.0/4.0)*Qi_iminus_n[i][k]-(1.0/8.0)*Qi_iminus_n[node[i].n_n[1]][k];  
						Qi_halfn_m[1][k] = (-1.0/8.0)*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]+(3.0/4.0)*Qi_iminus_n[node[i].n_n[3]][k]+(3.0/8.0)*Qi_iminus_n[i][k];
						Qi_halfn_m[2][k] = (3.0/8.0)*Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]-(5.0/4.0)*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]+(15.0/8.0)*Qi_iminus_n[node[i].n_n[3]][k];
						
						/** Qi(i-1/2)_plus**/				
						Qi_half_np[0][k] = (15.0/8.0)*Qi_iminus_n[i][k]-(5.0/4.0)*Qi_iminus_n[node[i].n_n[1]][k]+(3.0/8.0)*Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k];         
						Qi_half_np[1][k] = (3.0/8.0)*Qi_iminus_n[node[i].n_n[3]][k]+(3.0/4.0)*Qi_iminus_n[i][k]-(1.0/8.0)*Qi_iminus_n[node[i].n_n[1]][k];
						Qi_half_np[2][k] = (-1.0/8.0)*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]+(3.0/4.0)*Qi_iminus_n[node[i].n_n[3]][k]+(3.0/8.0)*Qi_iminus_n[i][k];
						
						/** Qj(i+1/2)_plus**/
						Qj_half_p[0][k] = (15.0/8.0)*Qj_iplus[node[i].n_n[0]][k]-(5.0/4.0)*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k]+(3.0/8.0)*Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k];           
						Qj_half_p[1][k] = (3.0/8.0)*Qj_iplus[i][k]+(3.0/4.0)*Qj_iplus[node[i].n_n[0]][k]-(1.0/8.0)*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k];
						Qj_half_p[2][k] = (-1.0/8.0)*Qj_iplus[node[i].n_n[2]][k]+(3.0/4.0)*Qj_iplus[i][k]+(3.0/8.0)*Qj_iplus[node[i].n_n[0]][k];
					  
						/** Qj(i+1/2)_minus**/
						Qj_half_m[0][k] = (3.0/8.0)*Qj_iplus[i][k]+(3.0/4.0)*Qj_iplus[node[i].n_n[0]][k]-(1.0/8.0)*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k];  
						Qj_half_m[1][k] = (-1.0/8.0)*Qj_iplus[node[i].n_n[2]][k]+(3.0/4.0)*Qj_iplus[i][k]+(3.0/8.0)*Qj_iplus[node[i].n_n[0]][k];
						Qj_half_m[2][k] = (3.0/8.0)*Qj_iplus[node[node[i].n_n[2]].n_n[2]][k]-(5.0/4.0)*Qj_iplus[node[i].n_n[2]][k]+(15.0/8.0)*Qj_iplus[i][k];
						
						/** Qj(i-1/2)_minus**/
						Qj_halfn_m[0][k] = (3.0/8.0)*Qj_iminus_n[node[i].n_n[2]][k]+(3.0/4.0)*Qj_iminus_n[i][k]-(1.0/8.0)*Qj_iminus_n[node[i].n_n[0]][k];  
						Qj_halfn_m[1][k] = (-1.0/8.0)*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]+(3.0/4.0)*Qj_iminus_n[node[i].n_n[2]][k]+(3.0/8.0)*Qj_iminus_n[i][k];
						Qj_halfn_m[2][k] = (3.0/8.0)*Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]-(5.0/4.0)*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]+(15.0/8.0)*Qj_iminus_n[node[i].n_n[2]][k];
						
						/** Qj(i-1/2)_plus**/				
						Qj_half_np[0][k] = (15.0/8.0)*Qj_iminus_n[i][k]-(5.0/4.0)*Qj_iminus_n[node[i].n_n[0]][k]+(3.0/8.0)*Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k];         
						Qj_half_np[1][k] = (3.0/8.0)*Qj_iminus_n[node[i].n_n[2]][k]+(3.0/4.0)*Qj_iminus_n[i][k]-(1.0/8.0)*Qj_iminus_n[node[i].n_n[0]][k];
						Qj_half_np[2][k] = (-1.0/8.0)*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]+(3.0/4.0)*Qj_iminus_n[node[i].n_n[2]][k]+(3.0/8.0)*Qj_iminus_n[i][k];
						
						/** Qk(i+1/2)_plus**/
						Qk_half_p[0][k] = (15.0/8.0)*Qk_iplus[node[i].n_n[4]][k]-(5.0/4.0)*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k]+(3.0/8.0)*Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k];           
						Qk_half_p[1][k] = (3.0/8.0)*Qk_iplus[i][k]+(3.0/4.0)*Qk_iplus[node[i].n_n[4]][k]-(1.0/8.0)*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k];
						Qk_half_p[2][k] = (-1.0/8.0)*Qk_iplus[node[i].n_n[5]][k]+(3.0/4.0)*Qk_iplus[i][k]+(3.0/8.0)*Qk_iplus[node[i].n_n[4]][k];
					  
						/** Qk(i+1/2)_minus**/
						Qk_half_m[0][k] = (3.0/8.0)*Qk_iplus[i][k]+(3.0/4.0)*Qk_iplus[node[i].n_n[4]][k]-(1.0/8.0)*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k];  
						Qk_half_m[1][k] = (-1.0/8.0)*Qk_iplus[node[i].n_n[5]][k]+(3.0/4.0)*Qk_iplus[i][k]+(3.0/8.0)*Qk_iplus[node[i].n_n[4]][k];
						Qk_half_m[2][k] = (3.0/8.0)*Qk_iplus[node[node[i].n_n[5]].n_n[5]][k]-(5.0/4.0)*Qk_iplus[node[i].n_n[5]][k]+(15.0/8.0)*Qk_iplus[i][k];
						
						/** Qk(i-1/2)_minus**/
						Qk_halfn_m[0][k] = (3.0/8.0)*Qk_iminus_n[node[i].n_n[5]][k]+(3.0/4.0)*Qk_iminus_n[i][k]-(1.0/8.0)*Qk_iminus_n[node[i].n_n[4]][k];  
						Qk_halfn_m[1][k] = (-1.0/8.0)*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]+(3.0/4.0)*Qk_iminus_n[node[i].n_n[5]][k]+(3.0/8.0)*Qk_iminus_n[i][k];
						Qk_halfn_m[2][k] = (3.0/8.0)*Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]-(5.0/4.0)*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]+(15.0/8.0)*Qk_iminus_n[node[i].n_n[5]][k];
						
						/** Qk(i-1/2)_plus**/				
						Qk_half_np[0][k] = (15.0/8.0)*Qk_iminus_n[i][k]-(5.0/4.0)*Qk_iminus_n[node[i].n_n[4]][k]+(3.0/8.0)*Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k];         
						Qk_half_np[1][k] = (3.0/8.0)*Qk_iminus_n[node[i].n_n[5]][k]+(3.0/4.0)*Qk_iminus_n[i][k]-(1.0/8.0)*Qk_iminus_n[node[i].n_n[4]][k];
						Qk_half_np[2][k] = (-1.0/8.0)*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]+(3.0/4.0)*Qk_iminus_n[node[i].n_n[5]][k]+(3.0/8.0)*Qk_iminus_n[i][k];

						/**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

						IS_Qim[2][k] = (13.0/12.0)*(pow(Qi_iplus[node[node[i].n_n[3]].n_n[3]][k]-2.0*Qi_iplus[node[i].n_n[3]][k]+Qi_iplus[i][k],2))+(1.0/4.0)*(pow(Qi_iplus[node[node[i].n_n[3]].n_n[3]][k]-4.0*Qi_iplus[node[i].n_n[3]][k]+3.0*Qi_iplus[i][k],2));	
						IS_Qim[1][k] = (13.0/12.0)*(pow(Qi_iplus[node[i].n_n[3]][k]-2.0*Qi_iplus[i][k]+Qi_iplus[node[i].n_n[1]][k],2))+(1.0/4.0)*(pow(Qi_iplus[node[i].n_n[3]][k]-Qi_iplus[node[i].n_n[1]][k],2));
						IS_Qim[0][k] = (13.0/12.0)*(pow(Qi_iplus[i][k]-2.0*Qi_iplus[node[i].n_n[1]][k]+Qi_iplus[node[node[i].n_n[1]].n_n[1]][k],2))+(1.0/4.0)*(pow(3.0*Qi_iplus[i][k]-4.0*Qi_iplus[node[i].n_n[1]][k]+Qi_iplus[node[node[i].n_n[1]].n_n[1]][k],2));
						
						IS_Qip[2][k] = (13.0/12.0)*(pow(Qi_iplus[node[i].n_n[3]][k]-2.0*Qi_iplus[i][k]+Qi_iplus[node[i].n_n[1]][k],2))+(1.0/4.0)*(pow(Qi_iplus[node[i].n_n[3]][k]-4.0*Qi_iplus[i][k]+3.0*Qi_iplus[node[i].n_n[1]][k],2));
						IS_Qip[1][k] = (13.0/12.0)*(pow(Qi_iplus[i][k]-2.0*Qi_iplus[node[i].n_n[1]][k]+Qi_iplus[node[node[i].n_n[1]].n_n[1]][k],2))+(1.0/4.0)*(pow(Qi_iplus[i][k]-Qi_iplus[node[node[i].n_n[1]].n_n[1]][k],2));
						IS_Qip[0][k] = (13.0/12.0)*(pow(Qi_iplus[node[i].n_n[1]][k]-2.0*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k]+Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k],2))+(1.0/4.0)*(pow(3.0*Qi_iplus[node[i].n_n[1]][k]-4.0*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k]+Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k],2));
						
						IS_Qinm[2][k] = (13.0/12.0)*(pow(Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]-2.0*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]+Qi_iminus_n[node[i].n_n[3]][k],2))+(1.0/4.0)*(pow(Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]-4.0*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]+3.0*Qi_iminus_n[node[i].n_n[3]][k],2));
						IS_Qinm[1][k] = (13.0/12.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]-2.0*Qi_iminus_n[node[i].n_n[3]][k]+Qi_iminus_n[i][k],2))+(1.0/4.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]-Qi_iminus_n[i][k],2));
						IS_Qinm[0][k] = (13.0/12.0)*(pow(Qi_iminus_n[node[i].n_n[3]][k]-2.0*Qi_iminus_n[i][k]+Qi_iminus_n[node[i].n_n[1]][k],2))+(1.0/4.0)*(pow(3.0*Qi_iminus_n[node[i].n_n[3]][k]-4.0*Qi_iminus_n[i][k]+Qi_iminus_n[node[i].n_n[1]][k],2));
						
						IS_Qinp[2][k] = (13.0/12.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]-2.0*Qi_iminus_n[node[i].n_n[3]][k]+Qi_iminus_n[i][k],2))+(1.0/4.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k]-4.0*Qi_iminus_n[node[i].n_n[3]][k]+3.0*Qi_iminus_n[i][k],2));
						IS_Qinp[1][k] = (13.0/12.0)*(pow(Qi_iminus_n[node[i].n_n[3]][k]-2.0*Qi_iminus_n[i][k]+Qi_iminus_n[node[i].n_n[1]][k],2))+(1.0/4.0)*(pow(Qi_iminus_n[node[i].n_n[3]][k]-Qi_iminus_n[node[i].n_n[1]][k],2));
						IS_Qinp[0][k] = (13.0/12.0)*(pow(Qi_iminus_n[i][k]-2.0*Qi_iminus_n[node[i].n_n[1]][k]+Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k],2))+(1.0/4.0)*(pow(3.0*Qi_iminus_n[i][k]-4.0*Qi_iminus_n[node[i].n_n[1]][k]+Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k],2));
						
						IS_Qjm[2][k] = (13.0/12.0)*(pow(Qj_iplus[node[node[i].n_n[2]].n_n[2]][k]-2.0*Qj_iplus[node[i].n_n[2]][k]+Qj_iplus[i][k],2))+(1.0/4.0)*(pow(Qj_iplus[node[node[i].n_n[2]].n_n[2]][k]-4.0*Qj_iplus[node[i].n_n[2]][k]+3.0*Qj_iplus[i][k],2));	
						IS_Qjm[1][k] = (13.0/12.0)*(pow(Qj_iplus[node[i].n_n[2]][k]-2.0*Qj_iplus[i][k]+Qj_iplus[node[i].n_n[0]][k],2))+(1.0/4.0)*(pow(Qj_iplus[node[i].n_n[2]][k]-Qj_iplus[node[i].n_n[0]][k],2));
						IS_Qjm[0][k] = (13.0/12.0)*(pow(Qj_iplus[i][k]-2.0*Qj_iplus[node[i].n_n[0]][k]+Qj_iplus[node[node[i].n_n[0]].n_n[0]][k],2))+(1.0/4.0)*(pow(3.0*Qj_iplus[i][k]-4.0*Qj_iplus[node[i].n_n[0]][k]+Qj_iplus[node[node[i].n_n[0]].n_n[0]][k],2));
						
						IS_Qjp[2][k] = (13.0/12.0)*(pow(Qj_iplus[node[i].n_n[2]][k]-2.0*Qj_iplus[i][k]+Qj_iplus[node[i].n_n[0]][k],2))+(1.0/4.0)*(pow(Qj_iplus[node[i].n_n[2]][k]-4.0*Qj_iplus[i][k]+3.0*Qj_iplus[node[i].n_n[0]][k],2));
						IS_Qjp[1][k] = (13.0/12.0)*(pow(Qj_iplus[i][k]-2.0*Qj_iplus[node[i].n_n[0]][k]+Qj_iplus[node[node[i].n_n[0]].n_n[0]][k],2))+(1.0/4.0)*(pow(Qj_iplus[i][k]-Qj_iplus[node[node[i].n_n[0]].n_n[0]][k],2));
						IS_Qjp[0][k] = (13.0/12.0)*(pow(Qj_iplus[node[i].n_n[0]][k]-2.0*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k]+Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k],2))+(1.0/4.0)*(pow(3.0*Qj_iplus[node[i].n_n[0]][k]-4.0*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k]+Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k],2));
						
						IS_Qjnm[2][k] = (13.0/12.0)*(pow(Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]-2.0*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]+Qj_iminus_n[node[i].n_n[2]][k],2))+(1.0/4.0)*(pow(Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]-4.0*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]+3.0*Qj_iminus_n[node[i].n_n[2]][k],2));
						IS_Qjnm[1][k] = (13.0/12.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]-2.0*Qj_iminus_n[node[i].n_n[2]][k]+Qj_iminus_n[i][k],2))+(1.0/4.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]-Qj_iminus_n[i][k],2));
						IS_Qjnm[0][k] = (13.0/12.0)*(pow(Qj_iminus_n[node[i].n_n[2]][k]-2.0*Qj_iminus_n[i][k]+Qj_iminus_n[node[i].n_n[0]][k],2))+(1.0/4.0)*(pow(3.0*Qj_iminus_n[node[i].n_n[2]][k]-4.0*Qj_iminus_n[i][k]+Qj_iminus_n[node[i].n_n[0]][k],2));
						
						IS_Qjnp[2][k] = (13.0/12.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]-2.0*Qj_iminus_n[node[i].n_n[2]][k]+Qj_iminus_n[i][k],2))+(1.0/4.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k]-4.0*Qj_iminus_n[node[i].n_n[2]][k]+3.0*Qj_iminus_n[i][k],2));
						IS_Qjnp[1][k] = (13.0/12.0)*(pow(Qj_iminus_n[node[i].n_n[2]][k]-2.0*Qj_iminus_n[i][k]+Qj_iminus_n[node[i].n_n[0]][k],2))+(1.0/4.0)*(pow(Qj_iminus_n[node[i].n_n[2]][k]-Qj_iminus_n[node[i].n_n[0]][k],2));
						IS_Qjnp[0][k] = (13.0/12.0)*(pow(Qj_iminus_n[i][k]-2.0*Qj_iminus_n[node[i].n_n[0]][k]+Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k],2))+(1.0/4.0)*(pow(3.0*Qj_iminus_n[i][k]-4.0*Qj_iminus_n[node[i].n_n[0]][k]+Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k],2));
						
						IS_Qkm[2][k] = (13.0/12.0)*(pow(Qk_iplus[node[node[i].n_n[5]].n_n[5]][k]-2.0*Qk_iplus[node[i].n_n[5]][k]+Qk_iplus[i][k],2))+(1.0/4.0)*(pow(Qk_iplus[node[node[i].n_n[5]].n_n[5]][k]-4.0*Qk_iplus[node[i].n_n[5]][k]+3.0*Qk_iplus[i][k],2));	
						IS_Qkm[1][k] = (13.0/12.0)*(pow(Qk_iplus[node[i].n_n[5]][k]-2.0*Qk_iplus[i][k]+Qk_iplus[node[i].n_n[4]][k],2))+(1.0/4.0)*(pow(Qk_iplus[node[i].n_n[5]][k]-Qk_iplus[node[i].n_n[4]][k],2));
						IS_Qkm[0][k] = (13.0/12.0)*(pow(Qk_iplus[i][k]-2.0*Qk_iplus[node[i].n_n[4]][k]+Qk_iplus[node[node[i].n_n[4]].n_n[4]][k],2))+(1.0/4.0)*(pow(3.0*Qk_iplus[i][k]-4.0*Qk_iplus[node[i].n_n[4]][k]+Qk_iplus[node[node[i].n_n[4]].n_n[4]][k],2));
						
						IS_Qkp[2][k] = (13.0/12.0)*(pow(Qk_iplus[node[i].n_n[5]][k]-2.0*Qk_iplus[i][k]+Qk_iplus[node[i].n_n[4]][k],2))+(1.0/4.0)*(pow(Qk_iplus[node[i].n_n[5]][k]-4.0*Qk_iplus[i][k]+3.0*Qk_iplus[node[i].n_n[4]][k],2));
						IS_Qkp[1][k] = (13.0/12.0)*(pow(Qk_iplus[i][k]-2.0*Qk_iplus[node[i].n_n[4]][k]+Qk_iplus[node[node[i].n_n[4]].n_n[4]][k],2))+(1.0/4.0)*(pow(Qk_iplus[i][k]-Qk_iplus[node[node[i].n_n[4]].n_n[4]][k],2));
						IS_Qkp[0][k] = (13.0/12.0)*(pow(Qk_iplus[node[i].n_n[4]][k]-2.0*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k]+Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k],2))+(1.0/4.0)*(pow(3.0*Qk_iplus[node[i].n_n[4]][k]-4.0*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k]+Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k],2));
						
						IS_Qknm[2][k] = (13.0/12.0)*(pow(Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]-2.0*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]+Qk_iminus_n[node[i].n_n[5]][k],2))+(1.0/4.0)*(pow(Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]-4.0*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]+3.0*Qk_iminus_n[node[i].n_n[5]][k],2));
						IS_Qknm[1][k] = (13.0/12.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]-2.0*Qk_iminus_n[node[i].n_n[5]][k]+Qk_iminus_n[i][k],2))+(1.0/4.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]-Qk_iminus_n[i][k],2));
						IS_Qknm[0][k] = (13.0/12.0)*(pow(Qk_iminus_n[node[i].n_n[5]][k]-2.0*Qk_iminus_n[i][k]+Qk_iminus_n[node[i].n_n[4]][k],2))+(1.0/4.0)*(pow(3.0*Qk_iminus_n[node[i].n_n[5]][k]-4.0*Qk_iminus_n[i][k]+Qk_iminus_n[node[i].n_n[4]][k],2));
						
						IS_Qknp[2][k] = (13.0/12.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]-2.0*Qk_iminus_n[node[i].n_n[5]][k]+Qk_iminus_n[i][k],2))+(1.0/4.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k]-4.0*Qk_iminus_n[node[i].n_n[5]][k]+3.0*Qk_iminus_n[i][k],2));
						IS_Qknp[1][k] = (13.0/12.0)*(pow(Qk_iminus_n[node[i].n_n[5]][k]-2.0*Qk_iminus_n[i][k]+Qk_iminus_n[node[i].n_n[4]][k],2))+(1.0/4.0)*(pow(Qk_iminus_n[node[i].n_n[5]][k]-Qk_iminus_n[node[i].n_n[4]][k],2));
						IS_Qknp[0][k] = (13.0/12.0)*(pow(Qk_iminus_n[i][k]-2.0*Qk_iminus_n[node[i].n_n[4]][k]+Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k],2))+(1.0/4.0)*(pow(3.0*Qk_iminus_n[i][k]-4.0*Qk_iminus_n[node[i].n_n[4]][k]+Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k],2));

						w_Qip[0][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qip[0][k]),2.0));
						w_Qip[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qip[1][k]),2.0));
						w_Qip[2][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qip[2][k]),2.0));
						
						w_Qim[0][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qim[0][k]),2.0));
						w_Qim[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qim[1][k]),2.0));
						w_Qim[2][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qim[2][k]),2.0));
						
						w_Qinp[0][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qinp[0][k]),2.0));
						w_Qinp[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qinp[1][k]),2.0));
						w_Qinp[2][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qinp[2][k]),2.0));
						
						w_Qinm[0][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qinm[0][k]),2.0));
						w_Qinm[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qinm[1][k]),2.0));
						w_Qinm[2][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qinm[2][k]),2.0));
						
						w_Qjp[0][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qjp[0][k]),2.0));
						w_Qjp[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qjp[1][k]),2.0));
						w_Qjp[2][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qjp[2][k]),2.0));
						
						w_Qjm[0][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qjm[0][k]),2.0));
						w_Qjm[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qjm[1][k]),2.0));
						w_Qjm[2][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qjm[2][k]),2.0));
						
						w_Qjnp[0][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qjnp[0][k]),2.0));
						w_Qjnp[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qjnp[1][k]),2.0));
						w_Qjnp[2][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qjnp[2][k]),2.0));
						
						w_Qjnm[0][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qjnm[0][k]),2.0));
						w_Qjnm[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qjnm[1][k]),2.0));
						w_Qjnm[2][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qjnm[2][k]),2.0));
		
						w_Qkp[0][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qkp[0][k]),2.0));
						w_Qkp[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qkp[1][k]),2.0));
						w_Qkp[2][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qkp[2][k]),2.0));
						
						w_Qkm[0][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qkm[0][k]),2.0));
						w_Qkm[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qkm[1][k]),2.0));
						w_Qkm[2][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qkm[2][k]),2.0));
						
						w_Qknp[0][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qknp[0][k]),2.0));
						w_Qknp[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qknp[1][k]),2.0));
						w_Qknp[2][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qknp[2][k]),2.0));
						
						w_Qknm[0][k] = (5.0/16.0)/(pow(Epsilon+(IS_Qknm[0][k]),2.0));
						w_Qknm[1][k] = (5.0/8.0)/(pow(Epsilon+(IS_Qknm[1][k]),2.0));
						w_Qknm[2][k] = (1.0/16.0)/(pow(Epsilon+(IS_Qknm[2][k]),2.0));
						
						W_Qip[0][k] = w_Qip[0][k]/(w_Qip[0][k]+w_Qip[1][k]+w_Qip[2][k]);
						W_Qip[1][k] = w_Qip[1][k]/(w_Qip[0][k]+w_Qip[1][k]+w_Qip[2][k]);
						W_Qip[2][k] = w_Qip[2][k]/(w_Qip[0][k]+w_Qip[1][k]+w_Qip[2][k]);
					
						W_Qinp[0][k] = w_Qinp[0][k]/(w_Qinp[0][k]+w_Qinp[1][k]+w_Qinp[2][k]);
						W_Qinp[1][k] = w_Qinp[1][k]/(w_Qinp[0][k]+w_Qinp[1][k]+w_Qinp[2][k]);
						W_Qinp[2][k] = w_Qinp[2][k]/(w_Qinp[0][k]+w_Qinp[1][k]+w_Qinp[2][k]);
					
						W_Qim[0][k] = w_Qim[0][k]/(w_Qim[0][k]+w_Qim[1][k]+w_Qim[2][k]);
						W_Qim[1][k] = w_Qim[1][k]/(w_Qim[0][k]+w_Qim[1][k]+w_Qim[2][k]);
						W_Qim[2][k] = w_Qim[2][k]/(w_Qim[0][k]+w_Qim[1][k]+w_Qim[2][k]);
					
						W_Qinm[0][k] = w_Qinm[0][k]/(w_Qinm[0][k]+w_Qinm[1][k]+w_Qinm[2][k]);
						W_Qinm[1][k] = w_Qinm[1][k]/(w_Qinm[0][k]+w_Qinm[1][k]+w_Qinm[2][k]);
						W_Qinm[2][k] = w_Qinm[2][k]/(w_Qinm[0][k]+w_Qinm[1][k]+w_Qinm[2][k]);
					
						W_Qjp[0][k] = w_Qjp[0][k]/(w_Qjp[0][k]+w_Qjp[1][k]+w_Qjp[2][k]);
						W_Qjp[1][k] = w_Qjp[1][k]/(w_Qjp[0][k]+w_Qjp[1][k]+w_Qjp[2][k]);
						W_Qjp[2][k] = w_Qjp[2][k]/(w_Qjp[0][k]+w_Qjp[1][k]+w_Qjp[2][k]);
					
						W_Qjnp[0][k] = w_Qjnp[0][k]/(w_Qjnp[0][k]+w_Qjnp[1][k]+w_Qjnp[2][k]);
						W_Qjnp[1][k] = w_Qjnp[1][k]/(w_Qjnp[0][k]+w_Qjnp[1][k]+w_Qjnp[2][k]);
						W_Qjnp[2][k] = w_Qjnp[2][k]/(w_Qjnp[0][k]+w_Qjnp[1][k]+w_Qjnp[2][k]);
					
						W_Qjm[0][k] = w_Qjm[0][k]/(w_Qjm[0][k]+w_Qjm[1][k]+w_Qjm[2][k]);
						W_Qjm[1][k] = w_Qjm[1][k]/(w_Qjm[0][k]+w_Qjm[1][k]+w_Qjm[2][k]);
						W_Qjm[2][k] = w_Qjm[2][k]/(w_Qjm[0][k]+w_Qjm[1][k]+w_Qjm[2][k]);
					
						W_Qjnm[0][k] = w_Qjnm[0][k]/(w_Qjnm[0][k]+w_Qjnm[1][k]+w_Qjnm[2][k]);
						W_Qjnm[1][k] = w_Qjnm[1][k]/(w_Qjnm[0][k]+w_Qjnm[1][k]+w_Qjnm[2][k]);
						W_Qjnm[2][k] = w_Qjnm[2][k]/(w_Qjnm[0][k]+w_Qjnm[1][k]+w_Qjnm[2][k]);
						
						W_Qkp[0][k] = w_Qkp[0][k]/(w_Qkp[0][k]+w_Qkp[1][k]+w_Qkp[2][k]);
						W_Qkp[1][k] = w_Qkp[1][k]/(w_Qkp[0][k]+w_Qkp[1][k]+w_Qkp[2][k]);
						W_Qkp[2][k] = w_Qkp[2][k]/(w_Qkp[0][k]+w_Qkp[1][k]+w_Qkp[2][k]);
					
						W_Qknp[0][k] = w_Qknp[0][k]/(w_Qknp[0][k]+w_Qknp[1][k]+w_Qknp[2][k]);
						W_Qknp[1][k] = w_Qknp[1][k]/(w_Qknp[0][k]+w_Qknp[1][k]+w_Qknp[2][k]);
						W_Qknp[2][k] = w_Qknp[2][k]/(w_Qknp[0][k]+w_Qknp[1][k]+w_Qknp[2][k]);
					
						W_Qkm[0][k] = w_Qkm[0][k]/(w_Qkm[0][k]+w_Qkm[1][k]+w_Qkm[2][k]);
						W_Qkm[1][k] = w_Qkm[1][k]/(w_Qkm[0][k]+w_Qkm[1][k]+w_Qkm[2][k]);
						W_Qkm[2][k] = w_Qkm[2][k]/(w_Qkm[0][k]+w_Qkm[1][k]+w_Qkm[2][k]);
					
						W_Qknm[0][k] = w_Qknm[0][k]/(w_Qknm[0][k]+w_Qknm[1][k]+w_Qknm[2][k]);
						W_Qknm[1][k] = w_Qknm[1][k]/(w_Qknm[0][k]+w_Qknm[1][k]+w_Qknm[2][k]);
						W_Qknm[2][k] = w_Qknm[2][k]/(w_Qknm[0][k]+w_Qknm[1][k]+w_Qknm[2][k]);
						
						Qi_iplus_half_pos_char[k] = W_Qip[0][k]*Qi_half_p[0][k]+W_Qip[1][k]*Qi_half_p[1][k]+W_Qip[2][k]*Qi_half_p[2][k];				
						Qi_iminus_half_pos_char[k] = W_Qinp[0][k]*Qi_half_np[0][k]+W_Qinp[1][k]*Qi_half_np[1][k]+W_Qinp[2][k]*Qi_half_np[2][k];	
					
						Qi_iplus_half_neg_char[k] = W_Qim[0][k]*Qi_half_m[0][k]+W_Qim[1][k]*Qi_half_m[1][k]+W_Qim[2][k]*Qi_half_m[2][k];				
						Qi_iminus_half_neg_char[k] = W_Qinm[0][k]*Qi_halfn_m[0][k]+W_Qinm[1][k]*Qi_halfn_m[1][k]+W_Qinm[2][k]*Qi_halfn_m[2][k];
					
						Qj_iplus_half_pos_char[k] = W_Qjp[0][k]*Qj_half_p[0][k]+W_Qjp[1][k]*Qj_half_p[1][k]+W_Qjp[2][k]*Qj_half_p[2][k];				
						Qj_iminus_half_pos_char[k] = W_Qjnp[0][k]*Qj_half_np[0][k]+W_Qjnp[1][k]*Qj_half_np[1][k]+W_Qjnp[2][k]*Qj_half_np[2][k];
					
						Qj_iplus_half_neg_char[k] = W_Qjm[0][k]*Qj_half_m[0][k]+W_Qjm[1][k]*Qj_half_m[1][k]+W_Qjm[2][k]*Qj_half_m[2][k];				
						Qj_iminus_half_neg_char[k] = W_Qjnm[0][k]*Qj_halfn_m[0][k]+W_Qjnm[1][k]*Qj_halfn_m[1][k]+W_Qjnm[2][k]*Qj_halfn_m[2][k];

						Qk_iplus_half_pos_char[k] = W_Qkp[0][k]*Qk_half_p[0][k]+W_Qkp[1][k]*Qk_half_p[1][k]+W_Qkp[2][k]*Qk_half_p[2][k];				
						Qk_iminus_half_pos_char[k] = W_Qknp[0][k]*Qk_half_np[0][k]+W_Qknp[1][k]*Qk_half_np[1][k]+W_Qknp[2][k]*Qk_half_np[2][k];
					
						Qk_iplus_half_neg_char[k] = W_Qkm[0][k]*Qk_half_m[0][k]+W_Qkm[1][k]*Qk_half_m[1][k]+W_Qkm[2][k]*Qk_half_m[2][k];				
						Qk_iminus_half_neg_char[k] = W_Qknm[0][k]*Qk_halfn_m[0][k]+W_Qknm[1][k]*Qk_halfn_m[1][k]+W_Qknm[2][k]*Qk_halfn_m[2][k];
						

						/***********************DELETE after checking******************************************************/
						/*********************************************************************************************************/
						/*********************************************************************************************************/
					/*
						Qi_iplus_half_pos[k] = W_Qip[0][k]*Qi_half_p[0][k]+W_Qip[1][k]*Qi_half_p[1][k]+W_Qip[2][k]*Qi_half_p[2][k];				
						Qi_iminus_half_pos[k] = W_Qinp[0][k]*Qi_half_np[0][k]+W_Qinp[1][k]*Qi_half_np[1][k]+W_Qinp[2][k]*Qi_half_np[2][k];	
					
						Qi_iplus_half_neg[k] = W_Qim[0][k]*Qi_half_m[0][k]+W_Qim[1][k]*Qi_half_m[1][k]+W_Qim[2][k]*Qi_half_m[2][k];				
						Qi_iminus_half_neg[k] = W_Qinm[0][k]*Qi_halfn_m[0][k]+W_Qinm[1][k]*Qi_halfn_m[1][k]+W_Qinm[2][k]*Qi_halfn_m[2][k];
					
						Qj_iplus_half_pos[k] = W_Qjp[0][k]*Qj_half_p[0][k]+W_Qjp[1][k]*Qj_half_p[1][k]+W_Qjp[2][k]*Qj_half_p[2][k];				
						Qj_iminus_half_pos[k] = W_Qjnp[0][k]*Qj_half_np[0][k]+W_Qjnp[1][k]*Qj_half_np[1][k]+W_Qjnp[2][k]*Qj_half_np[2][k];
					
						Qj_iplus_half_neg[k] = W_Qjm[0][k]*Qj_half_m[0][k]+W_Qjm[1][k]*Qj_half_m[1][k]+W_Qjm[2][k]*Qj_half_m[2][k];				
						Qj_iminus_half_neg[k] = W_Qjnm[0][k]*Qj_halfn_m[0][k]+W_Qjnm[1][k]*Qj_halfn_m[1][k]+W_Qjnm[2][k]*Qj_halfn_m[2][k];

						Qk_iplus_half_pos[k] = W_Qkp[0][k]*Qk_half_p[0][k]+W_Qkp[1][k]*Qk_half_p[1][k]+W_Qkp[2][k]*Qk_half_p[2][k];				
						Qk_iminus_half_pos[k] = W_Qknp[0][k]*Qk_half_np[0][k]+W_Qknp[1][k]*Qk_half_np[1][k]+W_Qknp[2][k]*Qk_half_np[2][k];
					
						Qk_iplus_half_neg[k] = W_Qkm[0][k]*Qk_half_m[0][k]+W_Qkm[1][k]*Qk_half_m[1][k]+W_Qkm[2][k]*Qk_half_m[2][k];				
						Qk_iminus_half_neg[k] = W_Qknm[0][k]*Qk_halfn_m[0][k]+W_Qknm[1][k]*Qk_halfn_m[1][k]+W_Qknm[2][k]*Qk_halfn_m[2][k];
						
						*/
						/*********************************************************************************************************/

					}
				
					m = 0;
					for(k=0;k<5;k++)
					{
						Qi_iplus_half_pos[k] = r_eigen_Qip[m][0]*(Qi_iplus_half_pos_char[0])+r_eigen_Qip[m][1]*(Qi_iplus_half_pos_char[1])+r_eigen_Qip[m][2]*(Qi_iplus_half_pos_char[2])+r_eigen_Qip[m][3]*(Qi_iplus_half_pos_char[3])+r_eigen_Qip[m][4]*(Qi_iplus_half_pos_char[4]);
						Qi_iplus_half_neg[k] = r_eigen_Qip[m][0]*(Qi_iplus_half_neg_char[0])+r_eigen_Qip[m][1]*(Qi_iplus_half_neg_char[1])+r_eigen_Qip[m][2]*(Qi_iplus_half_neg_char[2])+r_eigen_Qip[m][3]*(Qi_iplus_half_neg_char[3])+r_eigen_Qip[m][4]*(Qi_iplus_half_neg_char[4]);
						
						Qi_iminus_half_pos[k] = r_eigen_Qim[m][0]*(Qi_iminus_half_pos_char[0])+r_eigen_Qim[m][1]*(Qi_iminus_half_pos_char[1])+r_eigen_Qim[m][2]*(Qi_iminus_half_pos_char[2])+r_eigen_Qim[m][3]*(Qi_iminus_half_pos_char[3])+r_eigen_Qim[m][4]*(Qi_iminus_half_pos_char[4]);
						Qi_iminus_half_neg[k] = r_eigen_Qim[m][0]*(Qi_iminus_half_neg_char[0])+r_eigen_Qim[m][1]*(Qi_iminus_half_neg_char[1])+r_eigen_Qim[m][2]*(Qi_iminus_half_neg_char[2])+r_eigen_Qim[m][3]*(Qi_iminus_half_neg_char[3])+r_eigen_Qim[m][4]*(Qi_iminus_half_neg_char[4]);
							
						Qj_iplus_half_pos[k] = r_eigen_Qjp[m][0]*(Qj_iplus_half_pos_char[0])+r_eigen_Qjp[m][1]*(Qj_iplus_half_pos_char[1])+r_eigen_Qjp[m][2]*(Qj_iplus_half_pos_char[2])+r_eigen_Qjp[m][3]*(Qj_iplus_half_pos_char[3])+r_eigen_Qjp[m][4]*(Qj_iplus_half_pos_char[4]);
						Qj_iplus_half_neg[k] = r_eigen_Qjp[m][0]*(Qj_iplus_half_neg_char[0])+r_eigen_Qjp[m][1]*(Qj_iplus_half_neg_char[1])+r_eigen_Qjp[m][2]*(Qj_iplus_half_neg_char[2])+r_eigen_Qjp[m][3]*(Qj_iplus_half_neg_char[3])+r_eigen_Qjp[m][4]*(Qj_iplus_half_neg_char[4]);
						
						Qj_iminus_half_pos[k] = r_eigen_Qjm[m][0]*(Qj_iminus_half_pos_char[0])+r_eigen_Qjm[m][1]*(Qj_iminus_half_pos_char[1])+r_eigen_Qjm[m][2]*(Qj_iminus_half_pos_char[2])+r_eigen_Qjm[m][3]*(Qj_iminus_half_pos_char[3])+r_eigen_Qjm[m][4]*(Qj_iminus_half_pos_char[4]);
						Qj_iminus_half_neg[k] = r_eigen_Qjm[m][0]*(Qj_iminus_half_neg_char[0])+r_eigen_Qjm[m][1]*(Qj_iminus_half_neg_char[1])+r_eigen_Qjm[m][2]*(Qj_iminus_half_neg_char[2])+r_eigen_Qjm[m][3]*(Qj_iminus_half_neg_char[3])+r_eigen_Qjm[m][4]*(Qj_iminus_half_neg_char[4]);
						
						Qk_iplus_half_pos[k] = r_eigen_Qkp[m][0]*(Qk_iplus_half_pos_char[0])+r_eigen_Qkp[m][1]*(Qk_iplus_half_pos_char[1])+r_eigen_Qkp[m][2]*(Qk_iplus_half_pos_char[2])+r_eigen_Qkp[m][3]*(Qk_iplus_half_pos_char[3])+r_eigen_Qkp[m][4]*(Qk_iplus_half_pos_char[4]);
						Qk_iplus_half_neg[k] = r_eigen_Qkp[m][0]*(Qk_iplus_half_neg_char[0])+r_eigen_Qkp[m][1]*(Qk_iplus_half_neg_char[1])+r_eigen_Qkp[m][2]*(Qk_iplus_half_neg_char[2])+r_eigen_Qkp[m][3]*(Qk_iplus_half_neg_char[3])+r_eigen_Qkp[m][4]*(Qk_iplus_half_neg_char[4]);
						
						Qk_iminus_half_pos[k] = r_eigen_Qkm[m][0]*(Qk_iminus_half_pos_char[0])+r_eigen_Qkm[m][1]*(Qk_iminus_half_pos_char[1])+r_eigen_Qkm[m][2]*(Qk_iminus_half_pos_char[2])+r_eigen_Qkm[m][3]*(Qk_iminus_half_pos_char[3])+r_eigen_Qkm[m][4]*(Qk_iminus_half_pos_char[4]);
						Qk_iminus_half_neg[k] = r_eigen_Qkm[m][0]*(Qk_iminus_half_neg_char[0])+r_eigen_Qkm[m][1]*(Qk_iminus_half_neg_char[1])+r_eigen_Qkm[m][2]*(Qk_iminus_half_neg_char[2])+r_eigen_Qkm[m][3]*(Qk_iminus_half_neg_char[3])+r_eigen_Qkm[m][4]*(Qk_iminus_half_neg_char[4]);
					
						m++;
					}
				
					F_ip_pos[0] = Qi_iplus_half_pos[1];
					F_ip_pos[1] = (pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+0.4*(Qi_iplus_half_pos[4]-0.5*((pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0]))); 	
					F_ip_pos[2] = (Qi_iplus_half_pos[1]*Qi_iplus_half_pos[2])/Qi_iplus_half_pos[0];
					F_ip_pos[3] = (Qi_iplus_half_pos[1]*Qi_iplus_half_pos[3])/Qi_iplus_half_pos[0];
					F_ip_pos[4] = (1.4*Qi_iplus_half_pos[4]-0.2*((pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[1]/Qi_iplus_half_pos[0]);	
				
					F_ip_neg[0] = Qi_iplus_half_neg[1];
					F_ip_neg[1] = (pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+0.4*(Qi_iplus_half_neg[4]-0.5*((pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0]))); 	
					F_ip_neg[2] = (Qi_iplus_half_neg[1]*Qi_iplus_half_neg[2])/Qi_iplus_half_neg[0];
					F_ip_neg[3] = (Qi_iplus_half_neg[1]*Qi_iplus_half_neg[3])/Qi_iplus_half_neg[0];
					F_ip_neg[4] = (1.4*Qi_iplus_half_neg[4]-0.2*((pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[1]/Qi_iplus_half_neg[0]);	
				
					F_im_pos[0] = Qi_iminus_half_pos[1];
					F_im_pos[1] = (pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+0.4*(Qi_iminus_half_pos[4]-0.5*((pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0]))); 	
					F_im_pos[2] = (Qi_iminus_half_pos[1]*Qi_iminus_half_pos[2])/Qi_iminus_half_pos[0];
					F_im_pos[3] = (Qi_iminus_half_pos[1]*Qi_iminus_half_pos[3])/Qi_iminus_half_pos[0];
					F_im_pos[4] = (1.4*Qi_iminus_half_pos[4]-0.2*((pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[1]/Qi_iminus_half_pos[0]);	
				
					F_im_neg[0] = Qi_iminus_half_neg[1];
					F_im_neg[1] = (pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+0.4*(Qi_iminus_half_neg[4]-0.5*((pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0]))); 	
					F_im_neg[2] = (Qi_iminus_half_neg[1]*Qi_iminus_half_neg[2])/Qi_iminus_half_neg[0];
					F_im_neg[3] = (Qi_iminus_half_neg[1]*Qi_iminus_half_neg[3])/Qi_iminus_half_neg[0];
					F_im_neg[4] = (1.4*Qi_iminus_half_neg[4]-0.2*((pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[1]/Qi_iminus_half_neg[0]);	
				
					E_ip_pos[0] = Qi_iplus_half_pos[2];
					E_ip_pos[1] = (Qi_iplus_half_pos[1]*Qi_iplus_half_pos[2])/Qi_iplus_half_pos[0];
					E_ip_pos[2] = (pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+0.4*(Qi_iplus_half_pos[4]-0.5*((pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0]))); 	
					E_ip_pos[3] = (Qi_iplus_half_pos[2]*Qi_iplus_half_pos[3])/Qi_iplus_half_pos[0];
					E_ip_pos[4] = (1.4*Qi_iplus_half_pos[4]-0.2*((pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[2]/Qi_iplus_half_pos[0]);	
				
					E_ip_neg[0] = Qi_iplus_half_neg[2];
					E_ip_neg[1] = (Qi_iplus_half_neg[1]*Qi_iplus_half_neg[2])/Qi_iplus_half_neg[0];
					E_ip_neg[2] = (pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+0.4*(Qi_iplus_half_neg[4]-0.5*((pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0]))); 	
					E_ip_neg[3] = (Qi_iplus_half_neg[2]*Qi_iplus_half_neg[3])/Qi_iplus_half_neg[0];
					E_ip_neg[4] = (1.4*Qi_iplus_half_neg[4]-0.2*((pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[2]/Qi_iplus_half_neg[0]);	
						
					E_im_pos[0] = Qi_iminus_half_pos[2];
					E_im_pos[1] = (Qi_iminus_half_pos[1]*Qi_iminus_half_pos[2])/Qi_iminus_half_pos[0];
					E_im_pos[2] = (pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+0.4*(Qi_iminus_half_pos[4]-0.5*((pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0]))); 	
					E_im_pos[3] = (Qi_iminus_half_pos[2]*Qi_iminus_half_pos[3])/Qi_iminus_half_pos[0];
					E_im_pos[4] = (1.4*Qi_iminus_half_pos[4]-0.2*((pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[2]/Qi_iminus_half_pos[0]);	
						
					E_im_neg[0] = Qi_iminus_half_neg[2];
					E_im_neg[1] = (Qi_iminus_half_neg[1]*Qi_iminus_half_neg[2])/Qi_iminus_half_neg[0];
					E_im_neg[2] = (pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+0.4*(Qi_iminus_half_neg[4]-0.5*((pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0]))); 	
					E_im_neg[3] = (Qi_iminus_half_neg[2]*Qi_iminus_half_neg[3])/Qi_iminus_half_neg[0];
					E_im_neg[4] = (1.4*Qi_iminus_half_neg[4]-0.2*((pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[2]/Qi_iminus_half_neg[0]);
					
					G_ip_pos[0] = Qi_iplus_half_pos[3];
					G_ip_pos[1] = (Qi_iplus_half_pos[1]*Qi_iplus_half_pos[3])/Qi_iplus_half_pos[0];
					G_ip_pos[2] = (Qi_iplus_half_pos[2]*Qi_iplus_half_pos[3])/Qi_iplus_half_pos[0];
					G_ip_pos[3] = (pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0])+0.4*(Qi_iplus_half_pos[4]-0.5*((pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0]))); 	
					G_ip_pos[4] = (1.4*Qi_iplus_half_pos[4]-0.2*((pow(Qi_iplus_half_pos[1],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[2],2)/Qi_iplus_half_pos[0])+(pow(Qi_iplus_half_pos[3],2)/Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[3]/Qi_iplus_half_pos[0]);	
				
					G_ip_neg[0] = Qi_iplus_half_neg[3];
					G_ip_neg[1] = (Qi_iplus_half_neg[1]*Qi_iplus_half_neg[3])/Qi_iplus_half_neg[0];
					G_ip_neg[2] = (Qi_iplus_half_neg[2]*Qi_iplus_half_neg[3])/Qi_iplus_half_neg[0];
					G_ip_neg[3] = (pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0])+0.4*(Qi_iplus_half_neg[4]-0.5*((pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0]))); 	
					G_ip_neg[4] = (1.4*Qi_iplus_half_neg[4]-0.2*((pow(Qi_iplus_half_neg[1],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[2],2)/Qi_iplus_half_neg[0])+(pow(Qi_iplus_half_neg[3],2)/Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[3]/Qi_iplus_half_neg[0]);	
						
					G_im_pos[0] = Qi_iminus_half_pos[3];
					G_im_pos[1] = (Qi_iminus_half_pos[1]*Qi_iminus_half_pos[3])/Qi_iminus_half_pos[0];
					G_im_pos[2] = (Qi_iminus_half_pos[2]*Qi_iminus_half_pos[3])/Qi_iminus_half_pos[0];
					G_im_pos[3] = (pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0])+0.4*(Qi_iminus_half_pos[4]-0.5*((pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0]))); 	
					G_im_pos[4] = (1.4*Qi_iminus_half_pos[4]-0.2*((pow(Qi_iminus_half_pos[1],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[2],2)/Qi_iminus_half_pos[0])+(pow(Qi_iminus_half_pos[3],2)/Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[3]/Qi_iminus_half_pos[0]);	
						
					G_im_neg[0] = Qi_iminus_half_neg[3];
					G_im_neg[1] = (Qi_iminus_half_neg[1]*Qi_iminus_half_neg[3])/Qi_iminus_half_neg[0];
					G_im_neg[2] = (Qi_iminus_half_neg[2]*Qi_iminus_half_neg[3])/Qi_iminus_half_neg[0];
					G_im_neg[3] = (pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0])+0.4*(Qi_iminus_half_neg[4]-0.5*((pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0]))); 	
					G_im_neg[4] = (1.4*Qi_iminus_half_neg[4]-0.2*((pow(Qi_iminus_half_neg[1],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[2],2)/Qi_iminus_half_neg[0])+(pow(Qi_iminus_half_neg[3],2)/Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[3]/Qi_iminus_half_neg[0]);

					
					F_jp_pos[0] = Qj_iplus_half_pos[1];
					F_jp_pos[1] = (pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+0.4*(Qj_iplus_half_pos[4]-0.5*((pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0]))); 	
					F_jp_pos[2] = (Qj_iplus_half_pos[1]*Qj_iplus_half_pos[2])/Qj_iplus_half_pos[0];
					F_jp_pos[3] = (Qj_iplus_half_pos[1]*Qj_iplus_half_pos[3])/Qj_iplus_half_pos[0];
					F_jp_pos[4] = (1.4*Qj_iplus_half_pos[4]-0.2*((pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[1]/Qj_iplus_half_pos[0]);	
				
					F_jp_neg[0] = Qj_iplus_half_neg[1];
					F_jp_neg[1] = (pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+0.4*(Qj_iplus_half_neg[4]-0.5*((pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0]))); 	
					F_jp_neg[2] = (Qj_iplus_half_neg[1]*Qj_iplus_half_neg[2])/Qj_iplus_half_neg[0];
					F_jp_neg[3] = (Qj_iplus_half_neg[1]*Qj_iplus_half_neg[3])/Qj_iplus_half_neg[0];
					F_jp_neg[4] = (1.4*Qj_iplus_half_neg[4]-0.2*((pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[1]/Qj_iplus_half_neg[0]);	
				
					F_jm_pos[0] = Qj_iminus_half_pos[1];
					F_jm_pos[1] = (pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+0.4*(Qj_iminus_half_pos[4]-0.5*((pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0]))); 	
					F_jm_pos[2] = (Qj_iminus_half_pos[1]*Qj_iminus_half_pos[2])/Qj_iminus_half_pos[0];
					F_jm_pos[3] = (Qj_iminus_half_pos[1]*Qj_iminus_half_pos[3])/Qj_iminus_half_pos[0];
					F_jm_pos[4] = (1.4*Qj_iminus_half_pos[4]-0.2*((pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[1]/Qj_iminus_half_pos[0]);	
				
					F_jm_neg[0] = Qj_iminus_half_neg[1];
					F_jm_neg[1] = (pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+0.4*(Qj_iminus_half_neg[4]-0.5*((pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0]))); 	
					F_jm_neg[2] = (Qj_iminus_half_neg[1]*Qj_iminus_half_neg[2])/Qj_iminus_half_neg[0];
					F_jm_neg[3] = (Qj_iminus_half_neg[1]*Qj_iminus_half_neg[3])/Qj_iminus_half_neg[0];
					F_jm_neg[4] = (1.4*Qj_iminus_half_neg[4]-0.2*((pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[1]/Qj_iminus_half_neg[0]);	
				
					E_jp_pos[0] = Qj_iplus_half_pos[2];
					E_jp_pos[1] = (Qj_iplus_half_pos[1]*Qj_iplus_half_pos[2])/Qj_iplus_half_pos[0];
					E_jp_pos[2] = (pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+0.4*(Qj_iplus_half_pos[4]-0.5*((pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0]))); 	
					E_jp_pos[3] = (Qj_iplus_half_pos[2]*Qj_iplus_half_pos[3])/Qj_iplus_half_pos[0];
					E_jp_pos[4] = (1.4*Qj_iplus_half_pos[4]-0.2*((pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[2]/Qj_iplus_half_pos[0]);	
				
					E_jp_neg[0] = Qj_iplus_half_neg[2];
					E_jp_neg[1] = (Qj_iplus_half_neg[1]*Qj_iplus_half_neg[2])/Qj_iplus_half_neg[0];
					E_jp_neg[2] = (pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+0.4*(Qj_iplus_half_neg[4]-0.5*((pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0]))); 	
					E_jp_neg[3] = (Qj_iplus_half_neg[2]*Qj_iplus_half_neg[3])/Qj_iplus_half_neg[0];
					E_jp_neg[4] = (1.4*Qj_iplus_half_neg[4]-0.2*((pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[2]/Qj_iplus_half_neg[0]);	
						
					E_jm_pos[0] = Qj_iminus_half_pos[2];
					E_jm_pos[1] = (Qj_iminus_half_pos[1]*Qj_iminus_half_pos[2])/Qj_iminus_half_pos[0];
					E_jm_pos[2] = (pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+0.4*(Qj_iminus_half_pos[4]-0.5*((pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0]))); 	
					E_jm_pos[3] = (Qj_iminus_half_pos[2]*Qj_iminus_half_pos[3])/Qj_iminus_half_pos[0];
					E_jm_pos[4] = (1.4*Qj_iminus_half_pos[4]-0.2*((pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[2]/Qj_iminus_half_pos[0]);	
						
					E_jm_neg[0] = Qj_iminus_half_neg[2];
					E_jm_neg[1] = (Qj_iminus_half_neg[1]*Qj_iminus_half_neg[2])/Qj_iminus_half_neg[0];
					E_jm_neg[2] = (pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+0.4*(Qj_iminus_half_neg[4]-0.5*((pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0]))); 	
					E_jm_neg[3] = (Qj_iminus_half_neg[2]*Qj_iminus_half_neg[3])/Qj_iminus_half_neg[0];
					E_jm_neg[4] = (1.4*Qj_iminus_half_neg[4]-0.2*((pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[2]/Qj_iminus_half_neg[0]);
					
					G_jp_pos[0] = Qj_iplus_half_pos[3];
					G_jp_pos[1] = (Qj_iplus_half_pos[1]*Qj_iplus_half_pos[3])/Qj_iplus_half_pos[0];
					G_jp_pos[2] = (Qj_iplus_half_pos[2]*Qj_iplus_half_pos[3])/Qj_iplus_half_pos[0];
					G_jp_pos[3] = (pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0])+0.4*(Qj_iplus_half_pos[4]-0.5*((pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0]))); 	
					G_jp_pos[4] = (1.4*Qj_iplus_half_pos[4]-0.2*((pow(Qj_iplus_half_pos[1],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[2],2)/Qj_iplus_half_pos[0])+(pow(Qj_iplus_half_pos[3],2)/Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[3]/Qj_iplus_half_pos[0]);	
				
					G_jp_neg[0] = Qj_iplus_half_neg[3];
					G_jp_neg[1] = (Qj_iplus_half_neg[1]*Qj_iplus_half_neg[3])/Qj_iplus_half_neg[0];
					G_jp_neg[2] = (Qj_iplus_half_neg[2]*Qj_iplus_half_neg[3])/Qj_iplus_half_neg[0];
					G_jp_neg[3] = (pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0])+0.4*(Qj_iplus_half_neg[4]-0.5*((pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0]))); 	
					G_jp_neg[4] = (1.4*Qj_iplus_half_neg[4]-0.2*((pow(Qj_iplus_half_neg[1],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[2],2)/Qj_iplus_half_neg[0])+(pow(Qj_iplus_half_neg[3],2)/Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[3]/Qj_iplus_half_neg[0]);	
						
					G_jm_pos[0] = Qj_iminus_half_pos[3];
					G_jm_pos[1] = (Qj_iminus_half_pos[1]*Qj_iminus_half_pos[3])/Qj_iminus_half_pos[0];
					G_jm_pos[2] = (Qj_iminus_half_pos[2]*Qj_iminus_half_pos[3])/Qj_iminus_half_pos[0];
					G_jm_pos[3] = (pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0])+0.4*(Qj_iminus_half_pos[4]-0.5*((pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0]))); 	
					G_jm_pos[4] = (1.4*Qj_iminus_half_pos[4]-0.2*((pow(Qj_iminus_half_pos[1],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[2],2)/Qj_iminus_half_pos[0])+(pow(Qj_iminus_half_pos[3],2)/Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[3]/Qj_iminus_half_pos[0]);	
						
					G_jm_neg[0] = Qj_iminus_half_neg[3];
					G_jm_neg[1] = (Qj_iminus_half_neg[1]*Qj_iminus_half_neg[3])/Qj_iminus_half_neg[0];
					G_jm_neg[2] = (Qj_iminus_half_neg[2]*Qj_iminus_half_neg[3])/Qj_iminus_half_neg[0];
					G_jm_neg[3] = (pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0])+0.4*(Qj_iminus_half_neg[4]-0.5*((pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0]))); 	
					G_jm_neg[4] = (1.4*Qj_iminus_half_neg[4]-0.2*((pow(Qj_iminus_half_neg[1],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[2],2)/Qj_iminus_half_neg[0])+(pow(Qj_iminus_half_neg[3],2)/Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[3]/Qj_iminus_half_neg[0]);

					
					F_kp_pos[0] = Qk_iplus_half_pos[1];
					F_kp_pos[1] = (pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+0.4*(Qk_iplus_half_pos[4]-0.5*((pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0]))); 	
					F_kp_pos[2] = (Qk_iplus_half_pos[1]*Qk_iplus_half_pos[2])/Qk_iplus_half_pos[0];
					F_kp_pos[3] = (Qk_iplus_half_pos[1]*Qk_iplus_half_pos[3])/Qk_iplus_half_pos[0];
					F_kp_pos[4] = (1.4*Qk_iplus_half_pos[4]-0.2*((pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[1]/Qk_iplus_half_pos[0]);	
				
					F_kp_neg[0] = Qk_iplus_half_neg[1];
					F_kp_neg[1] = (pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+0.4*(Qk_iplus_half_neg[4]-0.5*((pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0]))); 	
					F_kp_neg[2] = (Qk_iplus_half_neg[1]*Qk_iplus_half_neg[2])/Qk_iplus_half_neg[0];
					F_kp_neg[3] = (Qk_iplus_half_neg[1]*Qk_iplus_half_neg[3])/Qk_iplus_half_neg[0];
					F_kp_neg[4] = (1.4*Qk_iplus_half_neg[4]-0.2*((pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[1]/Qk_iplus_half_neg[0]);	
				
					F_km_pos[0] = Qk_iminus_half_pos[1];
					F_km_pos[1] = (pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+0.4*(Qk_iminus_half_pos[4]-0.5*((pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0]))); 	
					F_km_pos[2] = (Qk_iminus_half_pos[1]*Qk_iminus_half_pos[2])/Qk_iminus_half_pos[0];
					F_km_pos[3] = (Qk_iminus_half_pos[1]*Qk_iminus_half_pos[3])/Qk_iminus_half_pos[0];
					F_km_pos[4] = (1.4*Qk_iminus_half_pos[4]-0.2*((pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[1]/Qk_iminus_half_pos[0]);	
				
					F_km_neg[0] = Qk_iminus_half_neg[1];
					F_km_neg[1] = (pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+0.4*(Qk_iminus_half_neg[4]-0.5*((pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0]))); 	
					F_km_neg[2] = (Qk_iminus_half_neg[1]*Qk_iminus_half_neg[2])/Qk_iminus_half_neg[0];
					F_km_neg[3] = (Qk_iminus_half_neg[1]*Qk_iminus_half_neg[3])/Qk_iminus_half_neg[0];
					F_km_neg[4] = (1.4*Qk_iminus_half_neg[4]-0.2*((pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[1]/Qk_iminus_half_neg[0]);	
				
					E_kp_pos[0] = Qk_iplus_half_pos[2];
					E_kp_pos[1] = (Qk_iplus_half_pos[1]*Qk_iplus_half_pos[2])/Qk_iplus_half_pos[0];
					E_kp_pos[2] = (pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+0.4*(Qk_iplus_half_pos[4]-0.5*((pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0]))); 	
					E_kp_pos[3] = (Qk_iplus_half_pos[2]*Qk_iplus_half_pos[3])/Qk_iplus_half_pos[0];
					E_kp_pos[4] = (1.4*Qk_iplus_half_pos[4]-0.2*((pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[2]/Qk_iplus_half_pos[0]);	
				
					E_kp_neg[0] = Qk_iplus_half_neg[2];
					E_kp_neg[1] = (Qk_iplus_half_neg[1]*Qk_iplus_half_neg[2])/Qk_iplus_half_neg[0];
					E_kp_neg[2] = (pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+0.4*(Qk_iplus_half_neg[4]-0.5*((pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0]))); 	
					E_kp_neg[3] = (Qk_iplus_half_neg[2]*Qk_iplus_half_neg[3])/Qk_iplus_half_neg[0];
					E_kp_neg[4] = (1.4*Qk_iplus_half_neg[4]-0.2*((pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[2]/Qk_iplus_half_neg[0]);	
						
					E_km_pos[0] = Qk_iminus_half_pos[2];
					E_km_pos[1] = (Qk_iminus_half_pos[1]*Qk_iminus_half_pos[2])/Qk_iminus_half_pos[0];
					E_km_pos[2] = (pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+0.4*(Qk_iminus_half_pos[4]-0.5*((pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0]))); 	
					E_km_pos[3] = (Qk_iminus_half_pos[2]*Qk_iminus_half_pos[3])/Qk_iminus_half_pos[0];
					E_km_pos[4] = (1.4*Qk_iminus_half_pos[4]-0.2*((pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[2]/Qk_iminus_half_pos[0]);	
						
					E_km_neg[0] = Qk_iminus_half_neg[2];
					E_km_neg[1] = (Qk_iminus_half_neg[1]*Qk_iminus_half_neg[2])/Qk_iminus_half_neg[0];
					E_km_neg[2] = (pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+0.4*(Qk_iminus_half_neg[4]-0.5*((pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0]))); 	
					E_km_neg[3] = (Qk_iminus_half_neg[2]*Qk_iminus_half_neg[3])/Qk_iminus_half_neg[0];
					E_km_neg[4] = (1.4*Qk_iminus_half_neg[4]-0.2*((pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[2]/Qk_iminus_half_neg[0]);
					
					G_kp_pos[0] = Qk_iplus_half_pos[3];
					G_kp_pos[1] = (Qk_iplus_half_pos[1]*Qk_iplus_half_pos[3])/Qk_iplus_half_pos[0];
					G_kp_pos[2] = (Qk_iplus_half_pos[2]*Qk_iplus_half_pos[3])/Qk_iplus_half_pos[0];
					G_kp_pos[3] = (pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0])+0.4*(Qk_iplus_half_pos[4]-0.5*((pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0]))); 	
					G_kp_pos[4] = (1.4*Qk_iplus_half_pos[4]-0.2*((pow(Qk_iplus_half_pos[1],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[2],2)/Qk_iplus_half_pos[0])+(pow(Qk_iplus_half_pos[3],2)/Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[3]/Qk_iplus_half_pos[0]);	
				
					G_kp_neg[0] = Qk_iplus_half_neg[3];
					G_kp_neg[1] = (Qk_iplus_half_neg[1]*Qk_iplus_half_neg[3])/Qk_iplus_half_neg[0];
					G_kp_neg[2] = (Qk_iplus_half_neg[2]*Qk_iplus_half_neg[3])/Qk_iplus_half_neg[0];
					G_kp_neg[3] = (pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0])+0.4*(Qk_iplus_half_neg[4]-0.5*((pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0]))); 	
					G_kp_neg[4] = (1.4*Qk_iplus_half_neg[4]-0.2*((pow(Qk_iplus_half_neg[1],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[2],2)/Qk_iplus_half_neg[0])+(pow(Qk_iplus_half_neg[3],2)/Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[3]/Qk_iplus_half_neg[0]);	
						
					G_km_pos[0] = Qk_iminus_half_pos[3];
					G_km_pos[1] = (Qk_iminus_half_pos[1]*Qk_iminus_half_pos[3])/Qk_iminus_half_pos[0];
					G_km_pos[2] = (Qk_iminus_half_pos[2]*Qk_iminus_half_pos[3])/Qk_iminus_half_pos[0];
					G_km_pos[3] = (pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0])+0.4*(Qk_iminus_half_pos[4]-0.5*((pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0]))); 	
					G_km_pos[4] = (1.4*Qk_iminus_half_pos[4]-0.2*((pow(Qk_iminus_half_pos[1],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[2],2)/Qk_iminus_half_pos[0])+(pow(Qk_iminus_half_pos[3],2)/Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[3]/Qk_iminus_half_pos[0]);	
						
					G_km_neg[0] = Qk_iminus_half_neg[3];
					G_km_neg[1] = (Qk_iminus_half_neg[1]*Qk_iminus_half_neg[3])/Qk_iminus_half_neg[0];
					G_km_neg[2] = (Qk_iminus_half_neg[2]*Qk_iminus_half_neg[3])/Qk_iminus_half_neg[0];
					G_km_neg[3] = (pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0])+0.4*(Qk_iminus_half_neg[4]-0.5*((pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0]))); 	
					G_km_neg[4] = (1.4*Qk_iminus_half_neg[4]-0.2*((pow(Qk_iminus_half_neg[1],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[2],2)/Qk_iminus_half_neg[0])+(pow(Qk_iminus_half_neg[3],2)/Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[3]/Qk_iminus_half_neg[0]);

					zeta_xip = metric[i].zeta_xip;
					zeta_yip = metric[i].zeta_yip;
					zeta_zip = metric[i].zeta_zip;
					eta_xjp = metric[i].eta_xjp;
					eta_yjp = metric[i].eta_yjp;	
					eta_zjp = metric[i].eta_zjp;	
					xi_xkp = metric[i].xi_xkp;
					xi_ykp = metric[i].xi_ykp;	
					xi_zkp = metric[i].xi_zkp;
				
					zeta_xim = metric[i].zeta_xim;
					zeta_yim = metric[i].zeta_yim;
					zeta_zim = metric[i].zeta_zim;
					eta_xjm = metric[i].eta_xjm;
					eta_yjm = metric[i].eta_yjm;	
					eta_zjm = metric[i].eta_zjm;	
					xi_xkm = metric[i].xi_xkm;
					xi_ykm = metric[i].xi_ykm;	
					xi_zkm = metric[i].xi_zkm;	
					
					eigen_Qip[0] =  (zeta_xip*(Qi_iplus_half_pos[1]/Qi_iplus_half_pos[0])+zeta_yip*(Qi_iplus_half_pos[2]/Qi_iplus_half_pos[0])+zeta_zip*(Qi_iplus_half_pos[3]/Qi_iplus_half_pos[0]))-sqrt(1.4*0.4*((Qi_iplus_half_pos[4]/Qi_iplus_half_pos[0])-((pow(Qi_iplus_half_pos[1],2)+pow(Qi_iplus_half_pos[2],2)+pow(Qi_iplus_half_pos[3],2))/(2.0*pow(Qi_iplus_half_pos[0],2))))*(zeta_xip*zeta_xip+zeta_yip*zeta_yip+zeta_zip*zeta_zip));
					eigen_Qip[1] =  (zeta_xip*(Qi_iplus_half_pos[1]/Qi_iplus_half_pos[0])+zeta_yip*(Qi_iplus_half_pos[2]/Qi_iplus_half_pos[0])+zeta_zip*(Qi_iplus_half_pos[3]/Qi_iplus_half_pos[0]));
					eigen_Qip[2] =  (zeta_xip*(Qi_iplus_half_pos[1]/Qi_iplus_half_pos[0])+zeta_yip*(Qi_iplus_half_pos[2]/Qi_iplus_half_pos[0])+zeta_zip*(Qi_iplus_half_pos[3]/Qi_iplus_half_pos[0]));
					eigen_Qip[3] =  (zeta_xip*(Qi_iplus_half_pos[1]/Qi_iplus_half_pos[0])+zeta_yip*(Qi_iplus_half_pos[2]/Qi_iplus_half_pos[0])+zeta_zip*(Qi_iplus_half_pos[3]/Qi_iplus_half_pos[0]));
					eigen_Qip[4] =  (zeta_xip*(Qi_iplus_half_pos[1]/Qi_iplus_half_pos[0])+zeta_yip*(Qi_iplus_half_pos[2]/Qi_iplus_half_pos[0])+zeta_zip*(Qi_iplus_half_pos[3]/Qi_iplus_half_pos[0]))+sqrt(1.4*0.4*((Qi_iplus_half_pos[4]/Qi_iplus_half_pos[0])-((pow(Qi_iplus_half_pos[1],2)+pow(Qi_iplus_half_pos[2],2)+pow(Qi_iplus_half_pos[3],2))/(2.0*pow(Qi_iplus_half_pos[0],2))))*(zeta_xip*zeta_xip+zeta_yip*zeta_yip+zeta_zip*zeta_zip));
					
					eigen_Qim[0] =  (zeta_xip*(Qi_iplus_half_neg[1]/Qi_iplus_half_neg[0])+zeta_yip*(Qi_iplus_half_neg[2]/Qi_iplus_half_neg[0])+zeta_zip*(Qi_iplus_half_neg[3]/Qi_iplus_half_neg[0]))-sqrt(1.4*0.4*((Qi_iplus_half_neg[4]/Qi_iplus_half_neg[0])-((pow(Qi_iplus_half_neg[1],2)+pow(Qi_iplus_half_neg[2],2)+pow(Qi_iplus_half_neg[3],2))/(2.0*pow(Qi_iplus_half_neg[0],2))))*(zeta_xip*zeta_xip+zeta_yip*zeta_yip+zeta_zip*zeta_zip));
					eigen_Qim[1] =  (zeta_xip*(Qi_iplus_half_neg[1]/Qi_iplus_half_neg[0])+zeta_yip*(Qi_iplus_half_neg[2]/Qi_iplus_half_neg[0])+zeta_zip*(Qi_iplus_half_neg[3]/Qi_iplus_half_neg[0]));
					eigen_Qim[2] =  (zeta_xip*(Qi_iplus_half_neg[1]/Qi_iplus_half_neg[0])+zeta_yip*(Qi_iplus_half_neg[2]/Qi_iplus_half_neg[0])+zeta_zip*(Qi_iplus_half_neg[3]/Qi_iplus_half_neg[0]));
					eigen_Qim[3] =  (zeta_xip*(Qi_iplus_half_neg[1]/Qi_iplus_half_neg[0])+zeta_yip*(Qi_iplus_half_neg[2]/Qi_iplus_half_neg[0])+zeta_zip*(Qi_iplus_half_neg[3]/Qi_iplus_half_neg[0]));
					eigen_Qim[4] =  (zeta_xip*(Qi_iplus_half_neg[1]/Qi_iplus_half_neg[0])+zeta_yip*(Qi_iplus_half_neg[2]/Qi_iplus_half_neg[0])+zeta_zip*(Qi_iplus_half_neg[3]/Qi_iplus_half_neg[0]))+sqrt(1.4*0.4*((Qi_iplus_half_neg[4]/Qi_iplus_half_neg[0])-((pow(Qi_iplus_half_neg[1],2)+pow(Qi_iplus_half_neg[2],2)+pow(Qi_iplus_half_neg[3],2))/(2.0*pow(Qi_iplus_half_neg[0],2))))*(zeta_xip*zeta_xip+zeta_yip*zeta_yip+zeta_zip*zeta_zip));
						
					eigen_Qinp[0] = (zeta_xim*(Qi_iminus_half_pos[1]/Qi_iminus_half_pos[0])+zeta_yim*(Qi_iminus_half_pos[2]/Qi_iminus_half_pos[0])+zeta_zim*(Qi_iminus_half_pos[3]/Qi_iminus_half_pos[0]))-sqrt(1.4*0.4*((Qi_iminus_half_pos[4]/Qi_iminus_half_pos[0])-((pow(Qi_iminus_half_pos[1],2)+pow(Qi_iminus_half_pos[2],2)+pow(Qi_iminus_half_pos[3],2))/(2.0*pow(Qi_iminus_half_pos[0],2))))*(zeta_xim*zeta_xim+zeta_yim*zeta_yim+zeta_zim*zeta_zim));
					eigen_Qinp[1] = (zeta_xim*(Qi_iminus_half_pos[1]/Qi_iminus_half_pos[0])+zeta_yim*(Qi_iminus_half_pos[2]/Qi_iminus_half_pos[0])+zeta_zim*(Qi_iminus_half_pos[3]/Qi_iminus_half_pos[0]));
					eigen_Qinp[2] = (zeta_xim*(Qi_iminus_half_pos[1]/Qi_iminus_half_pos[0])+zeta_yim*(Qi_iminus_half_pos[2]/Qi_iminus_half_pos[0])+zeta_zim*(Qi_iminus_half_pos[3]/Qi_iminus_half_pos[0]));
					eigen_Qinp[3] = (zeta_xim*(Qi_iminus_half_pos[1]/Qi_iminus_half_pos[0])+zeta_yim*(Qi_iminus_half_pos[2]/Qi_iminus_half_pos[0])+zeta_zim*(Qi_iminus_half_pos[3]/Qi_iminus_half_pos[0]));
					eigen_Qinp[4] = (zeta_xim*(Qi_iminus_half_pos[1]/Qi_iminus_half_pos[0])+zeta_yim*(Qi_iminus_half_pos[2]/Qi_iminus_half_pos[0])+zeta_zim*(Qi_iminus_half_pos[3]/Qi_iminus_half_pos[0]))+sqrt(1.4*0.4*((Qi_iminus_half_pos[4]/Qi_iminus_half_pos[0])-((pow(Qi_iminus_half_pos[1],2)+pow(Qi_iminus_half_pos[2],2)+pow(Qi_iminus_half_pos[3],2))/(2.0*pow(Qi_iminus_half_pos[0],2))))*(zeta_xim*zeta_xim+zeta_yim*zeta_yim+zeta_zim*zeta_zim));
					
					eigen_Qinm[0] = (zeta_xim*(Qi_iminus_half_neg[1]/Qi_iminus_half_neg[0])+zeta_yim*(Qi_iminus_half_neg[2]/Qi_iminus_half_neg[0])+zeta_zim*(Qi_iminus_half_neg[3]/Qi_iminus_half_neg[0]))-sqrt(1.4*0.4*((Qi_iminus_half_neg[4]/Qi_iminus_half_neg[0])-((pow(Qi_iminus_half_neg[1],2)+pow(Qi_iminus_half_neg[2],2)+pow(Qi_iminus_half_neg[3],2))/(2.0*pow(Qi_iminus_half_neg[0],2))))*(zeta_xim*zeta_xim+zeta_yim*zeta_yim+zeta_zim*zeta_zim));
					eigen_Qinm[1] = (zeta_xim*(Qi_iminus_half_neg[1]/Qi_iminus_half_neg[0])+zeta_yim*(Qi_iminus_half_neg[2]/Qi_iminus_half_neg[0])+zeta_zim*(Qi_iminus_half_neg[3]/Qi_iminus_half_neg[0]));
					eigen_Qinm[2] = (zeta_xim*(Qi_iminus_half_neg[1]/Qi_iminus_half_neg[0])+zeta_yim*(Qi_iminus_half_neg[2]/Qi_iminus_half_neg[0])+zeta_zim*(Qi_iminus_half_neg[3]/Qi_iminus_half_neg[0]));
					eigen_Qinm[3] = (zeta_xim*(Qi_iminus_half_neg[1]/Qi_iminus_half_neg[0])+zeta_yim*(Qi_iminus_half_neg[2]/Qi_iminus_half_neg[0])+zeta_zim*(Qi_iminus_half_neg[3]/Qi_iminus_half_neg[0]));
					eigen_Qinm[4] = (zeta_xim*(Qi_iminus_half_neg[1]/Qi_iminus_half_neg[0])+zeta_yim*(Qi_iminus_half_neg[2]/Qi_iminus_half_neg[0])+zeta_zim*(Qi_iminus_half_neg[3]/Qi_iminus_half_neg[0]))+sqrt(1.4*0.4*((Qi_iminus_half_neg[4]/Qi_iminus_half_neg[0])-((pow(Qi_iminus_half_neg[1],2)+pow(Qi_iminus_half_neg[2],2)+pow(Qi_iminus_half_neg[3],2))/(2.0*pow(Qi_iminus_half_neg[0],2))))*(zeta_xim*zeta_xim+zeta_yim*zeta_yim+zeta_zim*zeta_zim));
						
					eigen_Qjp[0] = (eta_xjp*(Qj_iplus_half_pos[1]/Qj_iplus_half_pos[0])+eta_yjp*(Qj_iplus_half_pos[2]/Qj_iplus_half_pos[0])+eta_zjp*(Qj_iplus_half_pos[3]/Qj_iplus_half_pos[0]))-sqrt(1.4*0.4*((Qj_iplus_half_pos[4]/Qj_iplus_half_pos[0])-((pow(Qj_iplus_half_pos[1],2)+pow(Qj_iplus_half_pos[2],2)+pow(Qj_iplus_half_pos[3],2))/(2.0*pow(Qj_iplus_half_pos[0],2))))*(eta_xjp*eta_xjp+eta_yjp*eta_yjp+eta_zjp*eta_zjp));
					eigen_Qjp[1] = (eta_xjp*(Qj_iplus_half_pos[1]/Qj_iplus_half_pos[0])+eta_yjp*(Qj_iplus_half_pos[2]/Qj_iplus_half_pos[0])+eta_zjp*(Qj_iplus_half_pos[3]/Qj_iplus_half_pos[0]));
					eigen_Qjp[2] = (eta_xjp*(Qj_iplus_half_pos[1]/Qj_iplus_half_pos[0])+eta_yjp*(Qj_iplus_half_pos[2]/Qj_iplus_half_pos[0])+eta_zjp*(Qj_iplus_half_pos[3]/Qj_iplus_half_pos[0]));
					eigen_Qjp[3] = (eta_xjp*(Qj_iplus_half_pos[1]/Qj_iplus_half_pos[0])+eta_yjp*(Qj_iplus_half_pos[2]/Qj_iplus_half_pos[0])+eta_zjp*(Qj_iplus_half_pos[3]/Qj_iplus_half_pos[0]));
					eigen_Qjp[4] = (eta_xjp*(Qj_iplus_half_pos[1]/Qj_iplus_half_pos[0])+eta_yjp*(Qj_iplus_half_pos[2]/Qj_iplus_half_pos[0])+eta_zjp*(Qj_iplus_half_pos[3]/Qj_iplus_half_pos[0]))+sqrt(1.4*0.4*((Qj_iplus_half_pos[4]/Qj_iplus_half_pos[0])-((pow(Qj_iplus_half_pos[1],2)+pow(Qj_iplus_half_pos[2],2)+pow(Qj_iplus_half_pos[3],2))/(2.0*pow(Qj_iplus_half_pos[0],2))))*(eta_xjp*eta_xjp+eta_yjp*eta_yjp+eta_zjp*eta_zjp));
					
					eigen_Qjm[0] = (eta_xjp*(Qj_iplus_half_neg[1]/Qj_iplus_half_neg[0])+eta_yjp*(Qj_iplus_half_neg[2]/Qj_iplus_half_neg[0])+eta_zjp*(Qj_iplus_half_neg[3]/Qj_iplus_half_neg[0]))-sqrt(1.4*0.4*((Qj_iplus_half_neg[4]/Qj_iplus_half_neg[0])-((pow(Qj_iplus_half_neg[1],2)+pow(Qj_iplus_half_neg[2],2)+pow(Qj_iplus_half_neg[3],2))/(2.0*pow(Qj_iplus_half_neg[0],2))))*(eta_xjp*eta_xjp+eta_yjp*eta_yjp+eta_zjp*eta_zjp));
					eigen_Qjm[1] = (eta_xjp*(Qj_iplus_half_neg[1]/Qj_iplus_half_neg[0])+eta_yjp*(Qj_iplus_half_neg[2]/Qj_iplus_half_neg[0])+eta_zjp*(Qj_iplus_half_neg[3]/Qj_iplus_half_neg[0]));
					eigen_Qjm[2] = (eta_xjp*(Qj_iplus_half_neg[1]/Qj_iplus_half_neg[0])+eta_yjp*(Qj_iplus_half_neg[2]/Qj_iplus_half_neg[0])+eta_zjp*(Qj_iplus_half_neg[3]/Qj_iplus_half_neg[0]));
					eigen_Qjm[3] = (eta_xjp*(Qj_iplus_half_neg[1]/Qj_iplus_half_neg[0])+eta_yjp*(Qj_iplus_half_neg[2]/Qj_iplus_half_neg[0])+eta_zjp*(Qj_iplus_half_neg[3]/Qj_iplus_half_neg[0]));
					eigen_Qjm[4] = (eta_xjp*(Qj_iplus_half_neg[1]/Qj_iplus_half_neg[0])+eta_yjp*(Qj_iplus_half_neg[2]/Qj_iplus_half_neg[0])+eta_zjp*(Qj_iplus_half_neg[3]/Qj_iplus_half_neg[0]))+sqrt(1.4*0.4*((Qj_iplus_half_neg[4]/Qj_iplus_half_neg[0])-((pow(Qj_iplus_half_neg[1],2)+pow(Qj_iplus_half_neg[2],2)+pow(Qj_iplus_half_neg[3],2))/(2.0*pow(Qj_iplus_half_neg[0],2))))*(eta_xjp*eta_xjp+eta_yjp*eta_yjp+eta_zjp*eta_zjp));
				
					eigen_Qjnp[0] = (eta_xjm*(Qj_iminus_half_pos[1]/Qj_iminus_half_pos[0])+eta_yjm*(Qj_iminus_half_pos[2]/Qj_iminus_half_pos[0])+eta_zjm*(Qj_iminus_half_pos[3]/Qj_iminus_half_pos[0]))-sqrt(1.4*0.4*((Qj_iminus_half_pos[4]/Qj_iminus_half_pos[0])-((pow(Qj_iminus_half_pos[1],2)+pow(Qj_iminus_half_pos[2],2)+pow(Qj_iminus_half_pos[3],2))/(2.0*pow(Qj_iminus_half_pos[0],2))))*(eta_xjm*eta_xjm+eta_yjm*eta_yjm+eta_zjm*eta_zjm));
					eigen_Qjnp[1] = (eta_xjm*(Qj_iminus_half_pos[1]/Qj_iminus_half_pos[0])+eta_yjm*(Qj_iminus_half_pos[2]/Qj_iminus_half_pos[0])+eta_zjm*(Qj_iminus_half_pos[3]/Qj_iminus_half_pos[0]));
					eigen_Qjnp[2] = (eta_xjm*(Qj_iminus_half_pos[1]/Qj_iminus_half_pos[0])+eta_yjm*(Qj_iminus_half_pos[2]/Qj_iminus_half_pos[0])+eta_zjm*(Qj_iminus_half_pos[3]/Qj_iminus_half_pos[0]));
					eigen_Qjnp[3] = (eta_xjm*(Qj_iminus_half_pos[1]/Qj_iminus_half_pos[0])+eta_yjm*(Qj_iminus_half_pos[2]/Qj_iminus_half_pos[0])+eta_zjm*(Qj_iminus_half_pos[3]/Qj_iminus_half_pos[0]));
					eigen_Qjnp[4] = (eta_xjm*(Qj_iminus_half_pos[1]/Qj_iminus_half_pos[0])+eta_yjm*(Qj_iminus_half_pos[2]/Qj_iminus_half_pos[0])+eta_zjm*(Qj_iminus_half_pos[3]/Qj_iminus_half_pos[0]))+sqrt(1.4*0.4*((Qj_iminus_half_pos[4]/Qj_iminus_half_pos[0])-((pow(Qj_iminus_half_pos[1],2)+pow(Qj_iminus_half_pos[2],2)+pow(Qj_iminus_half_pos[3],2))/(2.0*pow(Qj_iminus_half_pos[0],2))))*(eta_xjm*eta_xjm+eta_yjm*eta_yjm+eta_zjm*eta_zjm));
					
					eigen_Qjnm[0] = (eta_xjm*(Qj_iminus_half_neg[1]/Qj_iminus_half_neg[0])+eta_yjm*(Qj_iminus_half_neg[2]/Qj_iminus_half_neg[0])+eta_zjm*(Qj_iminus_half_neg[3]/Qj_iminus_half_neg[0]))-sqrt(1.4*0.4*((Qj_iminus_half_neg[4]/Qj_iminus_half_neg[0])-((pow(Qj_iminus_half_neg[1],2)+pow(Qj_iminus_half_neg[2],2)+pow(Qj_iminus_half_neg[3],2))/(2.0*pow(Qj_iminus_half_neg[0],2))))*(eta_xjm*eta_xjm+eta_yjm*eta_yjm+eta_zjm*eta_zjm));
					eigen_Qjnm[1] = (eta_xjm*(Qj_iminus_half_neg[1]/Qj_iminus_half_neg[0])+eta_yjm*(Qj_iminus_half_neg[2]/Qj_iminus_half_neg[0])+eta_zjm*(Qj_iminus_half_neg[3]/Qj_iminus_half_neg[0]));
					eigen_Qjnm[2] = (eta_xjm*(Qj_iminus_half_neg[1]/Qj_iminus_half_neg[0])+eta_yjm*(Qj_iminus_half_neg[2]/Qj_iminus_half_neg[0])+eta_zjm*(Qj_iminus_half_neg[3]/Qj_iminus_half_neg[0]));
					eigen_Qjnm[3] = (eta_xjm*(Qj_iminus_half_neg[1]/Qj_iminus_half_neg[0])+eta_yjm*(Qj_iminus_half_neg[2]/Qj_iminus_half_neg[0])+eta_zjm*(Qj_iminus_half_neg[3]/Qj_iminus_half_neg[0]));
					eigen_Qjnm[4] = (eta_xjm*(Qj_iminus_half_neg[1]/Qj_iminus_half_neg[0])+eta_yjm*(Qj_iminus_half_neg[2]/Qj_iminus_half_neg[0])+eta_zjm*(Qj_iminus_half_neg[3]/Qj_iminus_half_neg[0]))+sqrt(1.4*0.4*((Qj_iminus_half_neg[4]/Qj_iminus_half_neg[0])-((pow(Qj_iminus_half_neg[1],2)+pow(Qj_iminus_half_neg[2],2)+pow(Qj_iminus_half_neg[3],2))/(2.0*pow(Qj_iminus_half_neg[0],2))))*(eta_xjm*eta_xjm+eta_yjm*eta_yjm+eta_zjm*eta_zjm));
					
					eigen_Qkp[0] = (xi_xkp*(Qk_iplus_half_pos[1]/Qk_iplus_half_pos[0])+xi_ykp*(Qk_iplus_half_pos[2]/Qk_iplus_half_pos[0])+xi_zkp*(Qk_iplus_half_pos[3]/Qk_iplus_half_pos[0]))-sqrt(1.4*0.4*((Qk_iplus_half_pos[4]/Qk_iplus_half_pos[0])-((pow(Qk_iplus_half_pos[1],2)+pow(Qk_iplus_half_pos[2],2)+pow(Qk_iplus_half_pos[3],2))/(2.0*pow(Qk_iplus_half_pos[0],2))))*(xi_xkp*xi_xkp+xi_ykp*xi_ykp+xi_zkp*xi_zkp));
					eigen_Qkp[1] = (xi_xkp*(Qk_iplus_half_pos[1]/Qk_iplus_half_pos[0])+xi_ykp*(Qk_iplus_half_pos[2]/Qk_iplus_half_pos[0])+xi_zkp*(Qk_iplus_half_pos[3]/Qk_iplus_half_pos[0]));
					eigen_Qkp[2] = (xi_xkp*(Qk_iplus_half_pos[1]/Qk_iplus_half_pos[0])+xi_ykp*(Qk_iplus_half_pos[2]/Qk_iplus_half_pos[0])+xi_zkp*(Qk_iplus_half_pos[3]/Qk_iplus_half_pos[0]));
					eigen_Qkp[3] = (xi_xkp*(Qk_iplus_half_pos[1]/Qk_iplus_half_pos[0])+xi_ykp*(Qk_iplus_half_pos[2]/Qk_iplus_half_pos[0])+xi_zkp*(Qk_iplus_half_pos[3]/Qk_iplus_half_pos[0]));
					eigen_Qkp[4] = (xi_xkp*(Qk_iplus_half_pos[1]/Qk_iplus_half_pos[0])+xi_ykp*(Qk_iplus_half_pos[2]/Qk_iplus_half_pos[0])+xi_zkp*(Qk_iplus_half_pos[3]/Qk_iplus_half_pos[0]))+sqrt(1.4*0.4*((Qk_iplus_half_pos[4]/Qk_iplus_half_pos[0])-((pow(Qk_iplus_half_pos[1],2)+pow(Qk_iplus_half_pos[2],2)+pow(Qk_iplus_half_pos[3],2))/(2.0*pow(Qk_iplus_half_pos[0],2))))*(xi_xkp*xi_xkp+xi_ykp*xi_ykp+xi_zkp*xi_zkp));
					
					eigen_Qkm[0] = (xi_xkp*(Qk_iplus_half_neg[1]/Qk_iplus_half_neg[0])+xi_ykp*(Qk_iplus_half_neg[2]/Qk_iplus_half_neg[0])+xi_zkp*(Qk_iplus_half_neg[3]/Qk_iplus_half_neg[0]))-sqrt(1.4*0.4*((Qk_iplus_half_neg[4]/Qk_iplus_half_neg[0])-((pow(Qk_iplus_half_neg[1],2)+pow(Qk_iplus_half_neg[2],2)+pow(Qk_iplus_half_neg[3],2))/(2.0*pow(Qk_iplus_half_neg[0],2))))*(xi_xkp*xi_xkp+xi_ykp*xi_ykp+xi_zkp*xi_zkp));
					eigen_Qkm[1] = (xi_xkp*(Qk_iplus_half_neg[1]/Qk_iplus_half_neg[0])+xi_ykp*(Qk_iplus_half_neg[2]/Qk_iplus_half_neg[0])+xi_zkp*(Qk_iplus_half_neg[3]/Qk_iplus_half_neg[0]));
					eigen_Qkm[2] = (xi_xkp*(Qk_iplus_half_neg[1]/Qk_iplus_half_neg[0])+xi_ykp*(Qk_iplus_half_neg[2]/Qk_iplus_half_neg[0])+xi_zkp*(Qk_iplus_half_neg[3]/Qk_iplus_half_neg[0]));
					eigen_Qkm[3] = (xi_xkp*(Qk_iplus_half_neg[1]/Qk_iplus_half_neg[0])+xi_ykp*(Qk_iplus_half_neg[2]/Qk_iplus_half_neg[0])+xi_zkp*(Qk_iplus_half_neg[3]/Qk_iplus_half_neg[0]));
					eigen_Qkm[4] = (xi_xkp*(Qk_iplus_half_neg[1]/Qk_iplus_half_neg[0])+xi_ykp*(Qk_iplus_half_neg[2]/Qk_iplus_half_neg[0])+xi_zkp*(Qk_iplus_half_neg[3]/Qk_iplus_half_neg[0]))+sqrt(1.4*0.4*((Qk_iplus_half_neg[4]/Qk_iplus_half_neg[0])-((pow(Qk_iplus_half_neg[1],2)+pow(Qk_iplus_half_neg[2],2)+pow(Qk_iplus_half_neg[3],2))/(2.0*pow(Qk_iplus_half_neg[0],2))))*(xi_xkp*xi_xkp+xi_ykp*xi_ykp+xi_zkp*xi_zkp));
				
					eigen_Qknp[0] = (xi_xkm*(Qk_iminus_half_pos[1]/Qk_iminus_half_pos[0])+xi_ykm*(Qk_iminus_half_pos[2]/Qk_iminus_half_pos[0])+xi_zkm*(Qk_iminus_half_pos[3]/Qk_iminus_half_pos[0]))-sqrt(1.4*0.4*((Qk_iminus_half_pos[4]/Qk_iminus_half_pos[0])-((pow(Qk_iminus_half_pos[1],2)+pow(Qk_iminus_half_pos[2],2)+pow(Qk_iminus_half_pos[3],2))/(2.0*pow(Qk_iminus_half_pos[0],2))))*(xi_xkm*xi_xkm+xi_ykm*xi_ykm+xi_zkm*xi_zkm));
					eigen_Qknp[1] = (xi_xkm*(Qk_iminus_half_pos[1]/Qk_iminus_half_pos[0])+xi_ykm*(Qk_iminus_half_pos[2]/Qk_iminus_half_pos[0])+xi_zkm*(Qk_iminus_half_pos[3]/Qk_iminus_half_pos[0]));
					eigen_Qknp[2] = (xi_xkm*(Qk_iminus_half_pos[1]/Qk_iminus_half_pos[0])+xi_ykm*(Qk_iminus_half_pos[2]/Qk_iminus_half_pos[0])+xi_zkm*(Qk_iminus_half_pos[3]/Qk_iminus_half_pos[0]));
					eigen_Qknp[3] = (xi_xkm*(Qk_iminus_half_pos[1]/Qk_iminus_half_pos[0])+xi_ykm*(Qk_iminus_half_pos[2]/Qk_iminus_half_pos[0])+xi_zkm*(Qk_iminus_half_pos[3]/Qk_iminus_half_pos[0]));
					eigen_Qknp[4] = (xi_xkm*(Qk_iminus_half_pos[1]/Qk_iminus_half_pos[0])+xi_ykm*(Qk_iminus_half_pos[2]/Qk_iminus_half_pos[0])+xi_zkm*(Qk_iminus_half_pos[3]/Qk_iminus_half_pos[0]))+sqrt(1.4*0.4*((Qk_iminus_half_pos[4]/Qk_iminus_half_pos[0])-((pow(Qk_iminus_half_pos[1],2)+pow(Qk_iminus_half_pos[2],2)+pow(Qk_iminus_half_pos[3],2))/(2.0*pow(Qk_iminus_half_pos[0],2))))*(xi_xkm*xi_xkm+xi_ykm*xi_ykm+xi_zkm*xi_zkm));
					
					eigen_Qknm[0] = (xi_xkm*(Qk_iminus_half_neg[1]/Qk_iminus_half_neg[0])+xi_ykm*(Qk_iminus_half_neg[2]/Qk_iminus_half_neg[0])+xi_zkm*(Qk_iminus_half_neg[3]/Qk_iminus_half_neg[0]))-sqrt(1.4*0.4*((Qk_iminus_half_neg[4]/Qk_iminus_half_neg[0])-((pow(Qk_iminus_half_neg[1],2)+pow(Qk_iminus_half_neg[2],2)+pow(Qk_iminus_half_neg[3],2))/(2.0*pow(Qk_iminus_half_neg[0],2))))*(xi_xkm*xi_xkm+xi_ykm*xi_ykm+xi_zkm*xi_zkm));
					eigen_Qknm[1] = (xi_xkm*(Qk_iminus_half_neg[1]/Qk_iminus_half_neg[0])+xi_ykm*(Qk_iminus_half_neg[2]/Qk_iminus_half_neg[0])+xi_zkm*(Qk_iminus_half_neg[3]/Qk_iminus_half_neg[0]));
					eigen_Qknm[2] = (xi_xkm*(Qk_iminus_half_neg[1]/Qk_iminus_half_neg[0])+xi_ykm*(Qk_iminus_half_neg[2]/Qk_iminus_half_neg[0])+xi_zkm*(Qk_iminus_half_neg[3]/Qk_iminus_half_neg[0]));
					eigen_Qknm[3] = (xi_xkm*(Qk_iminus_half_neg[1]/Qk_iminus_half_neg[0])+xi_ykm*(Qk_iminus_half_neg[2]/Qk_iminus_half_neg[0])+xi_zkm*(Qk_iminus_half_neg[3]/Qk_iminus_half_neg[0]));
					eigen_Qknm[4] = (xi_xkm*(Qk_iminus_half_neg[1]/Qk_iminus_half_neg[0])+xi_ykm*(Qk_iminus_half_neg[2]/Qk_iminus_half_neg[0])+xi_zkm*(Qk_iminus_half_neg[3]/Qk_iminus_half_neg[0]))+sqrt(1.4*0.4*((Qk_iminus_half_neg[4]/Qk_iminus_half_neg[0])-((pow(Qk_iminus_half_neg[1],2)+pow(Qk_iminus_half_neg[2],2)+pow(Qk_iminus_half_neg[3],2))/(2.0*pow(Qk_iminus_half_neg[0],2))))*(xi_xkm*xi_xkm+xi_ykm*xi_ykm+xi_zkm*xi_zkm));

					for (k=0; k<5; k++)
					{
						F_ip_pos_comp[k] = deter[i].ip*(zeta_xip*F_ip_pos[k]+zeta_yip*E_ip_pos[k]+zeta_zip*G_ip_pos[k]);
						F_ip_neg_comp[k] = deter[i].ip*(zeta_xip*F_ip_neg[k]+zeta_yip*E_ip_neg[k]+zeta_zip*G_ip_neg[k]);
						
						F_im_pos_comp[k] = deter[i].im*(zeta_xim*F_im_pos[k]+zeta_yim*E_im_pos[k]+zeta_zim*G_im_pos[k]);
						F_im_neg_comp[k] = deter[i].im*(zeta_xim*F_im_neg[k]+zeta_yim*E_im_neg[k]+zeta_zim*G_im_neg[k]);
						
						E_ip_pos_comp[k] = deter[i].jp*(eta_xjp*F_jp_pos[k]+eta_yjp*E_jp_pos[k]+eta_zjp*G_jp_pos[k]);
						E_ip_neg_comp[k] = deter[i].jp*(eta_xjp*F_jp_neg[k]+eta_yjp*E_jp_neg[k]+eta_zjp*G_jp_neg[k]);
						
						E_im_pos_comp[k] = deter[i].jm*(eta_xjm*F_jm_pos[k]+eta_yjm*E_jm_pos[k]+eta_zjm*G_jm_pos[k]);
						E_im_neg_comp[k] = deter[i].jm*(eta_xjm*F_jm_neg[k]+eta_yjm*E_jm_neg[k]+eta_zjm*G_jm_neg[k]);
						
						G_ip_pos_comp[k] = deter[i].kp*(xi_xkp*F_kp_pos[k]+xi_ykp*E_kp_pos[k]+xi_zkp*G_kp_pos[k]);
						G_ip_neg_comp[k] = deter[i].kp*(xi_xkp*F_kp_neg[k]+xi_ykp*E_kp_neg[k]+xi_zkp*G_kp_neg[k]);
						
						G_im_pos_comp[k] = deter[i].km*(xi_xkm*F_km_pos[k]+xi_ykm*E_km_pos[k]+xi_zkm*G_km_pos[k]);
						G_im_neg_comp[k] = deter[i].km*(xi_xkm*F_km_neg[k]+xi_ykm*E_km_neg[k]+xi_zkm*G_km_neg[k]);
						
						Qi_iplus_half_pos[k] = deter[i].ip*Qi_iplus_half_pos[k];
						Qi_iplus_half_neg[k] = deter[i].ip*Qi_iplus_half_neg[k];
						
						Qi_iminus_half_pos[k] = deter[i].im*Qi_iminus_half_pos[k];
						Qi_iminus_half_neg[k] = deter[i].im*Qi_iminus_half_neg[k];
						
						Qj_iplus_half_pos[k] = deter[i].jp*Qj_iplus_half_pos[k];
						Qj_iplus_half_neg[k] = deter[i].jp*Qj_iplus_half_neg[k];
						
						Qj_iminus_half_pos[k] = deter[i].jm*Qj_iminus_half_pos[k];
						Qj_iminus_half_neg[k] = deter[i].jm*Qj_iminus_half_neg[k];
						
						Qk_iplus_half_pos[k] = deter[i].kp*Qk_iplus_half_pos[k];
						Qk_iplus_half_neg[k] = deter[i].kp*Qk_iplus_half_neg[k];
						
						Qk_iminus_half_pos[k] = deter[i].km*Qk_iminus_half_pos[k];
						Qk_iminus_half_neg[k] = deter[i].km*Qk_iminus_half_neg[k];
						
						/******************************************************************************************************************************************************/
						/******************************************************CENTRAL 4th ORDER VISCOUS TERMS*****************************************************************/
						
						if (node[i].loc > 0 && node[i].loc <= 52 &&  (node[node[i].n_n[1]].loc == 100 ))
						{
							dFv[k] = (1.0/6.0)*(11.0*Fv[i][k]-18.0*Fv[node[i].n_n[3]][k]+9.0*Fv[node[node[i].n_n[3]].n_n[3]][k]-2.0*Fv[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]);					
							dF[k] = (1.0/6.0)*(11.0*F[i][k]-18.0*F[node[i].n_n[3]][k]+9.0*F[node[node[i].n_n[3]].n_n[3]][k]-2.0*F[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]);
						}
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[3]].loc == 100 ))
						{
							dFv[k] = (1.0/6.0)*(-11.0*Fv[i][k]+18.0*Fv[node[i].n_n[1]][k]-9.0*Fv[node[node[i].n_n[1]].n_n[1]][k]+2.0*Fv[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
							dF[k] = (1.0/6.0)*(-11.0*F[i][k]+18.0*F[node[i].n_n[1]][k]-9.0*F[node[node[i].n_n[1]].n_n[1]][k]+2.0*F[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);					
						}
						if (node[i].loc == 0 && node[node[i].n_n[1]].loc > 0 && node[node[i].n_n[1]].loc <= 52 && (node[node[node[i].n_n[1]].n_n[1]].loc == 100 ))
						{
							dFv[k] = (1.0/6.0)*((2.0)*Fv[node[i].n_n[1]][k]+3.0*Fv[i][k]-6.0*Fv[node[i].n_n[3]][k]+Fv[node[node[i].n_n[3]].n_n[3]][k]);
							dF[k] = (1.0/6.0)*((2.0)*F[node[i].n_n[1]][k]+3.0*F[i][k]-6.0*F[node[i].n_n[3]][k]+F[node[node[i].n_n[3]].n_n[3]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[3]].loc > 0 && node[node[i].n_n[3]].loc <= 52 && (node[node[node[i].n_n[3]].n_n[3]].loc == 100 ))
						{
							dFv[k] = (1.0/6.0)*((-2.0)*Fv[node[i].n_n[3]][k]-3.0*Fv[i][k]+6.0*Fv[node[i].n_n[1]][k]-Fv[node[node[i].n_n[1]].n_n[1]][k]);
							dF[k] = (1.0/6.0)*((-2.0)*F[node[i].n_n[3]][k]-3.0*F[i][k]+6.0*F[node[i].n_n[1]][k]-F[node[node[i].n_n[1]].n_n[1]][k]);
						}
						if (node[i].loc == 0 && node[node[node[i].n_n[3]].n_n[3]].loc != 100 && node[node[i].n_n[3]].loc != 100 && node[node[node[i].n_n[1]].n_n[1]].loc != 100 && node[node[i].n_n[1]].loc != 100)
						{
							dFv[k] = (1.0/12.0)*(8.0*(Fv[node[i].n_n[1]][k]-Fv[node[i].n_n[3]][k])-(Fv[node[node[i].n_n[1]].n_n[1]][k]-Fv[node[node[i].n_n[3]].n_n[3]][k]));
							dF[k] = (1.0/12.0)*(8.0*(F[node[i].n_n[1]][k]-F[node[i].n_n[3]][k])-(F[node[node[i].n_n[1]].n_n[1]][k]-F[node[node[i].n_n[3]].n_n[3]][k]));
						}
						
						
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[0]].loc == 100 ))
						{
							dEv[k] = (1.0/6.0)*(11.0*Ev[i][k]-18.0*Ev[node[i].n_n[2]][k]+9.0*Ev[node[node[i].n_n[2]].n_n[2]][k]-2.0*Ev[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]);					
							dE[k] = (1.0/6.0)*(11.0*E[i][k]-18.0*E[node[i].n_n[2]][k]+9.0*E[node[node[i].n_n[2]].n_n[2]][k]-2.0*E[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]);					
						}
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[2]].loc == 100 ))
						{
							dEv[k] = (1.0/6.0)*(-11.0*Ev[i][k]+18.0*Ev[node[i].n_n[0]][k]-9.0*Ev[node[node[i].n_n[0]].n_n[0]][k]+2.0*Ev[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);					
							dE[k] = (1.0/6.0)*(-11.0*E[i][k]+18.0*E[node[i].n_n[0]][k]-9.0*E[node[node[i].n_n[0]].n_n[0]][k]+2.0*E[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);					
						}
						if (node[i].loc == 0 && node[node[i].n_n[0]].loc > 0 && node[node[i].n_n[0]].loc <= 52 && (node[node[node[i].n_n[0]].n_n[0]].loc == 100 ))
						{
							dEv[k] = (1.0/6.0)*(2.0*Ev[node[i].n_n[0]][k]+3.0*Ev[i][k]-6.0*Ev[node[i].n_n[2]][k]+Ev[node[node[i].n_n[2]].n_n[2]][k]);					
							dE[k] = (1.0/6.0)*(2.0*E[node[i].n_n[0]][k]+3.0*E[i][k]-6.0*E[node[i].n_n[2]][k]+E[node[node[i].n_n[2]].n_n[2]][k]);					
						}
						if (node[i].loc == 0 && node[node[i].n_n[2]].loc > 0 && node[node[i].n_n[2]].loc <= 52 && (node[node[node[i].n_n[2]].n_n[2]].loc == 100 ))
						{
							dEv[k] = (1.0/6.0)*(-2.0*Ev[node[i].n_n[2]][k]-3.0*Ev[i][k]+6.0*Ev[node[i].n_n[0]][k]-Ev[node[node[i].n_n[0]].n_n[0]][k]);
							dE[k] = (1.0/6.0)*(-2.0*E[node[i].n_n[2]][k]-3.0*E[i][k]+6.0*E[node[i].n_n[0]][k]-E[node[node[i].n_n[0]].n_n[0]][k]);
						}
						if (node[i].loc == 0 && node[node[node[i].n_n[2]].n_n[2]].loc != 100 && node[node[i].n_n[2]].loc != 100 && node[node[node[i].n_n[0]].n_n[0]].loc != 100 && node[node[i].n_n[0]].loc != 100)
						{
							dEv[k] = (1.0/12.0)*(8.0*(Ev[node[i].n_n[0]][k]-Ev[node[i].n_n[2]][k])-(Ev[node[node[i].n_n[0]].n_n[0]][k]-Ev[node[node[i].n_n[2]].n_n[2]][k])); 
							dE[k] = (1.0/12.0)*(8.0*(E[node[i].n_n[0]][k]-E[node[i].n_n[2]][k])-(E[node[node[i].n_n[0]].n_n[0]][k]-E[node[node[i].n_n[2]].n_n[2]][k])); 
						}


						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[4]].loc == 100 ))
						{
							dGv[k] = (1.0/6.0)*(11.0*Gv[i][k]-18.0*Gv[node[i].n_n[5]][k]+9.0*Gv[node[node[i].n_n[5]].n_n[5]][k]-2.0*Gv[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]);					
							dG[k] = (1.0/6.0)*(11.0*G[i][k]-18.0*G[node[i].n_n[5]][k]+9.0*G[node[node[i].n_n[5]].n_n[5]][k]-2.0*G[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]);			
						}
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[5]].loc == 100 ))
						{
							dGv[k] = (1.0/6.0)*(-11.0*Gv[i][k]+18.0*Gv[node[i].n_n[4]][k]-9.0*Gv[node[node[i].n_n[4]].n_n[4]][k]+2.0*Gv[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);					
							dG[k] = (1.0/6.0)*(-11.0*G[i][k]+18.0*G[node[i].n_n[4]][k]-9.0*G[node[node[i].n_n[4]].n_n[4]][k]+2.0*G[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);							
						}
						if (node[i].loc == 0 && node[node[i].n_n[4]].loc > 0 && node[node[i].n_n[4]].loc <= 52 && (node[node[node[i].n_n[4]].n_n[4]].loc == 100 ))
						{
							dGv[k] = (1.0/6.0)*(2.0*Gv[node[i].n_n[4]][k]+3.0*Gv[i][k]-6.0*Gv[node[i].n_n[5]][k]+Gv[node[node[i].n_n[5]].n_n[5]][k]);					
							dG[k] = (1.0/6.0)*(2.0*G[node[i].n_n[4]][k]+3.0*G[i][k]-6.0*G[node[i].n_n[5]][k]+G[node[node[i].n_n[5]].n_n[5]][k]);					
						}
						if (node[i].loc == 0 && node[node[i].n_n[5]].loc > 0 && node[node[i].n_n[5]].loc <= 52 && (node[node[node[i].n_n[5]].n_n[5]].loc == 100 ))
						{
							dGv[k] = (1.0/6.0)*(-2.0*Gv[node[i].n_n[5]][k]-3.0*Gv[i][k]+6.0*Gv[node[i].n_n[4]][k]-Gv[node[node[i].n_n[4]].n_n[4]][k]);
							dG[k] = (1.0/6.0)*(-2.0*G[node[i].n_n[5]][k]-3.0*G[i][k]+6.0*G[node[i].n_n[4]][k]-G[node[node[i].n_n[4]].n_n[4]][k]);
						}
						if (node[i].loc == 0 && node[node[node[i].n_n[5]].n_n[5]].loc != 100 && node[node[i].n_n[5]].loc != 100 && node[node[node[i].n_n[4]].n_n[4]].loc != 100 && node[node[i].n_n[4]].loc != 100)
						{
							dGv[k] = (1.0/12.0)*(8.0*(Gv[node[i].n_n[4]][k]-Gv[node[i].n_n[5]][k])-(Gv[node[node[i].n_n[4]].n_n[4]][k]-Gv[node[node[i].n_n[5]].n_n[5]][k])); 
							dG[k] = (1.0/12.0)*(8.0*(G[node[i].n_n[4]][k]-G[node[i].n_n[5]][k])-(G[node[node[i].n_n[4]].n_n[4]][k]-G[node[node[i].n_n[5]].n_n[5]][k])); 
						}
			//		}
					
			//		for (k=0; k<5; k++)
			//		{

						if (fabs(eigen_Qip[4]) >= fabs(eigen_Qim[4]))
						{
							alpha_u_ip[k] = fabs(eigen_Qip[4]);
						}
						else if (fabs(eigen_Qim[4]) >= fabs(eigen_Qip[4]))
						{
							alpha_u_ip[k] = fabs(eigen_Qim[4]);
						}
					
						if (fabs(eigen_Qjp[4]) >= fabs(eigen_Qjm[4]))
						{
							alpha_v_jp[k] = fabs(eigen_Qjp[4]);
						}
						else if (fabs(eigen_Qjm[4]) >= fabs(eigen_Qjp[4]))
						{
							alpha_v_jp[k] = fabs(eigen_Qjm[4]);
						}
						
						if (fabs(eigen_Qkp[4]) >= fabs(eigen_Qkm[4]))
						{
							alpha_w_kp[k] = fabs(eigen_Qkp[4]);
						}
						else if (fabs(eigen_Qkm[4]) >= fabs(eigen_Qkp[4]))
						{
							alpha_w_kp[k] = fabs(eigen_Qkm[4]);
						}
						
						if (fabs(eigen_Qinp[4]) >= fabs(eigen_Qinm[4]))
						{
							alpha_u_im[k] = fabs(eigen_Qinp[4]);
						}
						else if (fabs(eigen_Qinm[4]) >= fabs(eigen_Qinp[4]))
						{
							alpha_u_im[k] = fabs(eigen_Qinm[4]);
						}
						
						if (fabs(eigen_Qjnp[4]) >= fabs(eigen_Qjnm[4]))
						{
							alpha_v_jm[k] = fabs(eigen_Qjnp[4]);
						}
						else if (fabs(eigen_Qjnm[4]) >= fabs(eigen_Qjnp[4]))
						{
							alpha_v_jm[k] = fabs(eigen_Qjnm[4]);
						}
					
						if (fabs(eigen_Qknp[4]) >= fabs(eigen_Qknm[4]))
						{
							alpha_w_km[k] = fabs(eigen_Qknp[4]);
						}
						else if (fabs(eigen_Qknm[4]) >= fabs(eigen_Qknp[4]))
						{
							alpha_w_km[k] = fabs(eigen_Qknm[4]);
						}
						
						F_ip[k] = 0.5*(F_ip_pos_comp[k]+F_ip_neg_comp[k]-alpha_u_ip[k]*(Qi_iplus_half_pos[k]-Qi_iplus_half_neg[k]));	
						F_im[k] = 0.5*(F_im_pos_comp[k]+F_im_neg_comp[k]-alpha_u_im[k]*(Qi_iminus_half_pos[k]-Qi_iminus_half_neg[k]));	
						
						E_jp[k] = 0.5*(E_ip_pos_comp[k]+E_ip_neg_comp[k]-alpha_v_jp[k]*(Qj_iplus_half_pos[k]-Qj_iplus_half_neg[k]));		
						E_jm[k] = 0.5*(E_im_pos_comp[k]+E_im_neg_comp[k]-alpha_v_jm[k]*(Qj_iminus_half_pos[k]-Qj_iminus_half_neg[k]));	
						
						G_kp[k] = 0.5*(G_ip_pos_comp[k]+G_ip_neg_comp[k]-alpha_w_kp[k]*(Qk_iplus_half_pos[k]-Qk_iplus_half_neg[k]));		
						G_km[k] = 0.5*(G_im_pos_comp[k]+G_im_neg_comp[k]-alpha_w_km[k]*(Qk_iminus_half_pos[k]-Qk_iminus_half_neg[k]));	
				
						h_F_ip[k] = F_ip[k];
						h_F_im[k] = F_im[k];
						
						h_E_ip[k] = E_jp[k];
						h_E_im[k] = E_jm[k];
						
						h_G_ip[k] = G_kp[k];
						h_G_im[k] = G_km[k];
						
						d2F_d2z_ip[k] = (1.0/48.0)*(-5.0*F[node[node[i].n_n[3]].n_n[3]][k]+39.0*F[node[i].n_n[3]][k]-34.0*F[i][k]-34.0*F[node[i].n_n[1]][k]+39.0*F[node[node[i].n_n[1]].n_n[1]][k]-5.0*F[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
						d2F_d2z_im[k] = (1.0/48.0)*(-5.0*F[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]+39.0*F[node[node[i].n_n[3]].n_n[3]][k]-34.0*F[node[i].n_n[3]][k]-34.0*F[i][k]+39.0*F[node[i].n_n[1]][k]-5.0*F[node[node[i].n_n[1]].n_n[1]][k]);
						
						d2E_d2e_ip[k] = (1.0/48.0)*(-5.0*E[node[node[i].n_n[2]].n_n[2]][k]+39.0*E[node[i].n_n[2]][k]-34.0*E[i][k]-34.0*E[node[i].n_n[0]][k]+39.0*E[node[node[i].n_n[0]].n_n[0]][k]-5.0*E[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);
						d2E_d2e_im[k] = (1.0/48.0)*(-5.0*E[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]+39.0*E[node[node[i].n_n[2]].n_n[2]][k]-34.0*E[node[i].n_n[2]][k]-34.0*E[i][k]+39.0*E[node[i].n_n[0]][k]-5.0*E[node[node[i].n_n[0]].n_n[0]][k]);
						
						d2G_d2x_ip[k] = (1.0/48.0)*(-5.0*G[node[node[i].n_n[5]].n_n[5]][k]+39.0*G[node[i].n_n[5]][k]-34.0*G[i][k]-34.0*G[node[i].n_n[4]][k]+39.0*G[node[node[i].n_n[4]].n_n[4]][k]-5.0*G[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);
						d2G_d2x_im[k] = (1.0/48.0)*(-5.0*G[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]+39.0*G[node[node[i].n_n[5]].n_n[5]][k]-34.0*G[node[i].n_n[5]][k]-34.0*G[i][k]+39.0*G[node[i].n_n[4]][k]-5.0*G[node[node[i].n_n[4]].n_n[4]][k]);
						
						d4F_d4z_ip[k] = (1.0/2.0)*(F[node[node[i].n_n[3]].n_n[3]][k]-3.0*F[node[i].n_n[3]][k]+2.0*F[i][k]+2.0*F[node[i].n_n[1]][k]-3.0*F[node[node[i].n_n[1]].n_n[1]][k]+F[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
						d4F_d4z_im[k] = (1.0/2.0)*(F[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]-3.0*F[node[node[i].n_n[3]].n_n[3]][k]+2.0*F[node[i].n_n[3]][k]+2.0*F[i][k]-3.0*F[node[i].n_n[1]][k]+F[node[node[i].n_n[1]].n_n[1]][k]);
						
						d4E_d4e_ip[k] = (1.0/2.0)*(E[node[node[i].n_n[2]].n_n[2]][k]-3.0*E[node[i].n_n[2]][k]+2.0*E[i][k]+2.0*E[node[i].n_n[0]][k]-3.0*E[node[node[i].n_n[0]].n_n[0]][k]+E[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);
						d4E_d4e_im[k] = (1.0/2.0)*(E[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]-3.0*E[node[node[i].n_n[2]].n_n[2]][k]+2.0*E[node[i].n_n[2]][k]+2.0*E[i][k]-3.0*E[node[i].n_n[0]][k]+E[node[node[i].n_n[0]].n_n[0]][k]);
						
						d4G_d4x_ip[k] = (1.0/2.0)*(G[node[node[i].n_n[5]].n_n[5]][k]-3.0*G[node[i].n_n[5]][k]+2.0*G[i][k]+2.0*G[node[i].n_n[4]][k]-3.0*G[node[node[i].n_n[4]].n_n[4]][k]+G[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);
						d4G_d4x_im[k] = (1.0/2.0)*(G[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]-3.0*G[node[node[i].n_n[5]].n_n[5]][k]+2.0*G[node[i].n_n[5]][k]+2.0*G[i][k]-3.0*G[node[i].n_n[4]][k]+G[node[node[i].n_n[4]].n_n[4]][k]);
					
						F_ip_h[k] = h_F_ip[k]-(1.0/24.0)*d2F_d2z_ip[k]+(7.0/5760.0)*d4F_d4z_ip[k];
						F_im_h[k] = h_F_im[k]-(1.0/24.0)*d2F_d2z_im[k]+(7.0/5760.0)*d4F_d4z_im[k];
					
						E_ip_h[k] = h_E_ip[k]-(1.0/24.0)*d2E_d2e_ip[k]+(7.0/5760.0)*d4E_d4e_ip[k];
						E_im_h[k] = h_E_im[k]-(1.0/24.0)*d2E_d2e_im[k]+(7.0/5760.0)*d4E_d4e_im[k]; 
						
						G_ip_h[k] = h_G_ip[k]-(1.0/24.0)*d2G_d2x_ip[k]+(7.0/5760.0)*d4G_d4x_ip[k];
						G_im_h[k] = h_G_im[k]-(1.0/24.0)*d2G_d2x_im[k]+(7.0/5760.0)*d4G_d4x_im[k]; 
					
						dF_W[k] = F_ip_h[k]-F_im_h[k];
						dE_W[k] = E_ip_h[k]-E_im_h[k];
						dG_W[k] = G_ip_h[k]-G_im_h[k];
						
						dF_C[k] = dF[k];
						dE_C[k] = dE[k];
						dG_C[k] = dG[k];
						
						dF[k] = DUCROS[i]*dF_W[k]+(1.0-DUCROS[i])*dF_C[k];
						dE[k] = DUCROS[i]*dE_W[k]+(1.0-DUCROS[i])*dE_C[k];
						dG[k] = DUCROS[i]*dG_W[k]+(1.0-DUCROS[i])*dG_C[k];	
					
						dF[k] = F_ip_h[k]-F_im_h[k];
						dE[k] = E_ip_h[k]-E_im_h[k];
						dG[k] = G_ip_h[k]-G_im_h[k];
					
						L[i][j][k] = (-1.0)*((1.0/del_zeta)*dF[k]+(1.0/del_eta)*dE[k]+(1.0/del_eta)*dG[k]+(1.0/del_zeta)*dFv[k]+(1.0/del_eta)*dEv[k]+(1.0/del_eta)*dGv[k]);		
						/***********************3rd order RK Scheme***************************************************************************/
						if(j==0)
						{
							U[i][gar][k] = U[i][0][k]+del_t*L[i][0][k];
						}
						else if (j==1)
						{
							U[i][gar][k] = (3.0/4.0)*U[i][0][k]+(1.0/4.0)*U[i][1][k]+0.25*del_t*L[i][1][k];					
						}
						else if (j==2)
						{
							U[i][gar][k] = (U[i][0][k]/3.0)+((2.0/3.0)*U[i][2][k])+(2.0/3.0)*(del_t*L[i][2][k]);
							lm =0;
						}
						final_U[k] = (1.0/det[i])*U[i][gar][k];
					}
					
				//	fprintf(fnode,"%e %e %e %e %e %e %d %d\n", dFv[4], dEv[4], dGv[4], dF[4], dE[4], dG[4], node[i].global, node[i].loc);
					
					rho[lm][i] = final_U[0];
					u[lm][i] = final_U[1]/final_U[0];
					v[lm][i] = final_U[2]/final_U[0];
					w[lm][i] = final_U[3]/final_U[0];
					e[lm][i] = (final_U[4]/final_U[0])-0.5*(pow(u[lm][i],2.0)+pow(v[lm][i],2.0)+pow(w[lm][i],2.0));
					p[lm][i] = e[lm][i]*rho[lm][i]*0.4;
					t[lm][i] = (1.4*Mach*Mach)*(p[lm][i]/rho[lm][i]);
					a[lm][i] = sqrt(1.4*p[lm][i]/rho[lm][i]);
					
					if(j == 2)
					{	
						div[i][1].u = u[0][i];
						div[i][1].v = v[0][i];
						div[i][1].e = e[0][i];	
					
						if( max_div_u < fabs(div[i][1].u-div[i][0].u))
						{
							max_div_u = fabs(div[i][1].u-div[i][0].u);
							u_loc = i;
						}
						if( max_div_v < fabs(div[i][1].v-div[i][0].v))
						{
							max_div_v = fabs(div[i][1].v-div[i][0].v);
							v_loc = i;
						}
						if( max_div_e < fabs(div[i][1].e-div[i][0].e))
						{
							max_div_e = fabs(div[i][1].e-div[i][0].e);
							e_loc = i;
						}
					
						div[i][0].u = div[i][1].u;
						div[i][0].v = div[i][1].v;
						div[i][0].e = div[i][1].e;		
					}
				}
			}
			
		//	fclose(fnode);
			
			position = 0;
			for(i=0; i<neigh_pro; i++)
			{
				MPI_Pack_size(recv_c[c[i]]*9,MPI_DOUBLE,MPI_COMM_WORLD,&memsize1);
				buffer2 = malloc(memsize1);                                          /***********carefull with buffer1 and buffer2******************/
				MPI_Pack_size(proc_node[c[i]]*9,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
				buffer = malloc(memsize);
				position = 0;
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&u[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&v[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&w[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&p[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&t[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&rho[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&e[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for(k=0; k<recv_c[c[i]];k++)
				{
					MPI_Pack(&mu[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
			
				MPI_Sendrecv(buffer2,memsize1,MPI_PACKED,c[i],c[i],buffer,memsize,MPI_PACKED,c[i],myrank, MPI_COMM_WORLD, &status);
				position = 0;
				
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&u[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&v[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&w[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&p[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&t[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&rho[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&e[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				for(k=0; k<proc_node[c[i]];k++)
				{
					MPI_Unpack(buffer,memsize,&position,&mu[lm][recv_b[i][k]],1,MPI_DOUBLE,MPI_COMM_WORLD);
				}
				free(buffer);
				buffer = NULL;
				free(buffer2);
				buffer2 = NULL;
				k =0;
			}
	
			if (iter % 200000 == 0)
			{
				back_pressure = back_pressure+0.01;
			}
			initial =1;
			intialise(u, v, w, rho, p, t, mu, e, g_node, inlet_node, outlet_node, wall_node, boundary_node, CD, wal_node, initial, jacobian, metric, lm, a, \
			nodes, NUMNP, out_node, inl_node, node, bou_node, all_bou_node, all_boundary_nodes, div, sd_node, sd_wal_node, sd_bou_node, sd_out_node, sd_inl_node, sd_wall_node, \
			sd_boundary_node, sd_outlet_node, sd_inlet_node, singular);		
			
			for (k=1; k<=no_of_tip_send; k++)
			{
				u[lm][tip[k].node] = 0.0;
				v[lm][tip[k].node] = 0.0;
				w[lm][tip[k].node] = 0.0;
				p[lm][tip[k].node] = p[lm][node[tip[k].node].n_n[3]];
				t[lm][tip[k].node] = t[lm][node[tip[k].node].n_n[3]];
				rho[lm][tip[k].node] = (1.4*Mach*Mach)*(p[lm][tip[k].node]/t[lm][tip[k].node]);
				e[lm][tip[k].node] = p[lm][tip[k].node]/(0.4*rho[lm][tip[k].node]);
			}
			
			for (k=1; k<=no_of_tip_recv; k++)
			{	
				u[lm][tip_recv[k].node] = 0.0;
				v[lm][tip_recv[k].node] = 0.0;
				w[lm][tip_recv[k].node] = 0.0;
				p[lm][tip_recv[k].node] = p[lm][node[tip_recv[k].node].n_n[3]];
				t[lm][tip_recv[k].node] = t[lm][node[tip_recv[k].node].n_n[3]];
				rho[lm][tip_recv[k].node] = (1.4*Mach*Mach)*(p[lm][tip_recv[k].node]/t[lm][tip_recv[k].node]);
				e[lm][tip_recv[k].node] = p[lm][tip_recv[k].node]/(0.4*rho[lm][tip_recv[k].node]);
			}
		}
		
		end = time(NULL);			

		position = 0;
		MPI_Pack_size(3, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		MPI_Pack(&max_div_u, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Pack(&max_div_v, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Pack(&max_div_e, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		free(buffer);
		buffer = NULL;
		
		if ((iter%100)==0)
		{
			position =0;
			MPI_Pack_size(sd_node*10, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);
			MPI_Pack(&u[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&v[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&w[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&rho[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&p[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&t[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&e[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&mu[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&DUCROS[1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		
			MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);			
			free(buffer);
			buffer = NULL;	
		}
	}
}








