#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern double del_zeta;
extern double del_eta;
extern double Mach;
extern double Free_t;
extern float Free_rho;
extern float Free_mu;
extern float Free_pres;
extern float Free_kine;
extern float Free_sound;
extern float univ_gas;
extern float char_length;
//extern float gamma;
extern double Pr;
extern double del_t;
extern double Epsilon;
extern int restart_num;
extern double Reyl;
extern int iterations;
extern char *file_ext;
extern double back_pressure;

typedef struct 
{
	double x, y, z;
	int e[8];	
	int ID;
	int n_n[7];
//	int location;
//	int corner;
//	int weno;
	int corner_ID;
	int val;
	int loc;
	int proc;
	int local;
	int global;
	int req;
} MNODE;

typedef struct 
{
	int e[10];
} TMP_NODE;

typedef struct
{
	int n_n[7];
} singul;
	
typedef struct
{
	double u, v, e;
} RESIDUAL;	 

typedef struct
{
	int ID;
} COR_NODE;

/*
typedef struct 
{
	double zeta_x, zeta_y, zeta_z, eta_x, eta_y, eta_z, xi_x, xi_y, xi_z;
	double zeta_xip, zeta_yip, zeta_zip, eta_xip, eta_yip, eta_zip, xi_xip, xi_yip, xi_zip, zeta_xim, zeta_yim, zeta_zim, eta_xim, eta_yim, eta_zim, xi_xim, xi_yim, xi_zim;
	double zeta_xjp, zeta_yjp, zeta_zjp, eta_xjp, eta_yjp, eta_zjp, xi_xjp, xi_yjp, xi_zjp, zeta_xjm, zeta_yjm, zeta_zjm, eta_xjm, eta_yjm, eta_zjm, xi_xjm, xi_yjm, xi_zjm;
	double zeta_xkp, zeta_ykp, zeta_zkp, eta_xkp, eta_ykp, eta_zkp, xi_xkp, xi_ykp, xi_zkp, zeta_xkm, zeta_ykm, zeta_zkm, eta_xkm, eta_ykm, eta_zkm, xi_xkm, xi_ykm, xi_zkm;
		
	//double x_zeta_ip, y_zeta_ip, z_zeta_ip, x_eta_ip, y_eta_ip, z_eta_ip, x_xi_ip, y_xi_ip, z_xi_ip, x_zeta_im, y_zeta_im, z_zeta_im, x_eta_im, y_eta_im, z_eta_im, x_xi_im, y_xi_im, z_xi_im;
	//double x_zeta_jp, y_zeta_jp, z_zeta_jp, x_eta_jp, y_eta_jp, z_eta_jp, x_xi_jp, y_xi_jp, z_xi_jp, x_zeta_jm, y_zeta_jm, z_zeta_jm, x_eta_jm, y_eta_jm, z_eta_jm, x_xi_jm, y_xi_jm, z_xi_jm;
	//double x_zeta_kp, y_zeta_kp, z_zeta_kp, x_eta_kp, y_eta_kp, z_eta_kp, x_xi_kp, y_xi_kp, z_xi_kp, x_zeta_km, y_zeta_km, z_zeta_km, x_eta_km, y_eta_km, z_eta_km, x_xi_km, y_xi_km, z_xi_km;
	double ip, im, jp, jm, kp, km;
	double x_zeta, x_eta, x_xi, y_zeta, y_eta, y_xi, z_zeta, z_eta, z_xi;
} JACOB;

*/

typedef struct 
{
	double zeta_x, zeta_y, zeta_z, eta_x, eta_y, eta_z, xi_x, xi_y, xi_z;
	double zeta_xip, zeta_yip, zeta_zip, eta_xip, eta_yip, eta_zip, xi_xip, xi_yip, xi_zip, zeta_xim, zeta_yim, zeta_zim, eta_xim, eta_yim, eta_zim, xi_xim, xi_yim, xi_zim;
	double zeta_xjp, zeta_yjp, zeta_zjp, eta_xjp, eta_yjp, eta_zjp, xi_xjp, xi_yjp, xi_zjp, zeta_xjm, zeta_yjm, zeta_zjm, eta_xjm, eta_yjm, eta_zjm, xi_xjm, xi_yjm, xi_zjm;
	double zeta_xkp, zeta_ykp, zeta_zkp, eta_xkp, eta_ykp, eta_zkp, xi_xkp, xi_ykp, xi_zkp, zeta_xkm, zeta_ykm, zeta_zkm, eta_xkm, eta_ykm, eta_zkm, xi_xkm, xi_ykm, xi_zkm;
} TRANSFORMED;

typedef struct 
{
	double x_zeta, x_eta, x_xi, y_zeta, y_eta, y_xi, z_zeta, z_eta, z_xi;
} JACOB;

typedef struct 
{
	double ip;
} Qip;

typedef struct 
{
	double ip, im, jp, jm, kp, km;
} DETERM;

typedef struct
{
	double zeta_1,zeta_2,zeta_3; 
	double eta_1,eta_2,eta_3;
	double xi_1,xi_2,xi_3;
} JAC;

typedef struct
{
	double zeta_1,zeta_2,zeta_3,zeta_4,zeta_5,zeta_6; 
	double eta_1,eta_2,eta_3,eta_4,eta_5,eta_6;
	double xi_1,xi_2,xi_3,xi_4,xi_5,xi_6;
} TEMP_JAC;

typedef struct
{
	int connect[9], element_neighbour[9];
	int type;
	int location;
	int corner;
} ELEM;
	
typedef struct
	{
	int check;
	} BOUND;

typedef struct 
{
	int proc;
	int n_proc;
	int global;
} SUB_DOM ;	

typedef struct 
{
	int proc;
	int node;
	int i;
	int j;
	int k;
} DOM_TIP ;	


void metric_term(ELEM *CD, MNODE *node, int NELEM, JACOB *jacobian, TRANSFORMED *metric, int cor, int NUMNP, double *det, int *outlet_node, int out_node, int *boundary_node, int bou_node, int *inlet_node, int inl_node, int *wall_node, int wal_node, int *all_boundary_nodes, int all_bou_node, JAC *jac, TEMP_JAC *jaco);

void writepltfile(int NELEM, int NUMNP, MNODE *node, ELEM *CD, double **u, double **v, double **w, double **rho, double **p, double **t, double **e, double **mu, int iteration, double *DUCROS);

void restart_file(int NELEM, int NUMNP, MNODE *node, ELEM *CD, double **u, double **v, double **w, double **rho, double **p, double **t, double **e, double **mu, int iteration, int restart_num);

void writepltfile_jacob(int NELEM, int NUMNP, MNODE *node, ELEM *CD, double **u, double **v, double **rho, double **p, double **t, JACOB *jacobian, TRANSFORMED *metric, double *det);

void viscousflux_variables(ELEM *CD, double **u, double **v, double **w, double **t, int j, int NUMNP, MNODE *node, double **mu, double Reyl, \
double *tauzz, double *tauze, double *tauee, double *tauxx, double *tauzx, double *tauex, double *qz, double *qe, double *qx, \
TRANSFORMED *metric, int g_node, int i, int *wall_node, double ***U, singul *singular, double *TAU_SGS_XY, double *TAU_SGS_YZ, double *TAU_SGS_XZ,\
double *H_SGS_X, double *H_SGS_Y, double *H_SGS_Z, double *D_SGS_X, double *D_SGS_Y, double *D_SGS_Z, double *TAU_SGS_XX, double *TAU_SGS_YY, double *TAU_SGS_ZZ, double *det, double *DUCROS);

void roe_average(double **u, double **v, double **w, double **rho, double **p, double **t, double **mu, double **e, int g_node, ELEM *CD, JACOB *jacobian, TRANSFORMED *metric, \
int j, double **a, MNODE *node, double *roe_rho_ip, double *roe_u_ip, double *roe_v_ip, double *roe_w_ip, double *roe_h_ip, double *roe_rho_im, double *roe_u_im, double *roe_v_im, \
double *roe_w_im, double *roe_h_im, double *roe_rho_jp, double *roe_u_jp, double *roe_v_jp, double *roe_w_jp, double *roe_h_jp, double *roe_rho_jm, \
double *roe_u_jm, double *roe_v_jm, double *roe_w_jm, double *roe_h_jm, double *roe_rho_kp, double *roe_u_kp, double *roe_v_kp, double *roe_w_kp, double *roe_h_kp, \
double *roe_rho_km, double *roe_u_km, double *roe_v_km, double *roe_w_km, double *roe_h_km, int all_bou_node, int *all_boundary_nodes, RESIDUAL **div, int sd_node, int i); 

void intialise(double **u, double **v, double **w, double **rho, double **p, double **t, double **mu, double **e, int g_node, int *inlet_node, int *outlet_node, int *wall_node, \
int *boundary_node, ELEM *CD, int wal_node, int initial, JACOB *jacobian, TRANSFORMED *metric,int j, double **a, int nodes, int NUMNP, int out_node, int inl_node, MNODE *node, \
int bou_node, int all_bou_node, int *all_boundary_nodes, RESIDUAL **div, int sd_node, int sd_wal_node,int sd_bou_node, int sd_out_node, int sd_inl_node, int *sd_wall_node, \
int *sd_boundary_node, int *sd_outlet_node, int *sd_inlet_node, singul *singular);

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
double *H_SGS_X, double *H_SGS_Y, double *H_SGS_Z, double *D_SGS_X, double *D_SGS_Y, double *D_SGS_Z, double *TAU_SGS_XX, double *TAU_SGS_YY, double *TAU_SGS_ZZ, double *DUCROS);

void weno_solver_0(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus);

void weno_solver_1(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus);

void weno_solver_2(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus);

void weno_solver_3(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus);

void weno_solver_4(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus);

void weno_solver_5(ELEM *CD, double *roe_u, double *roe_v, double *roe_w, double *roe_h, double *roe_rho, double **r_eigen_f, double **l_eigen_f, \
Qip **Q, int j, int hk, int g_elem, double **u, double **v, double **w, double **e, double **p, double **rho, double ***U, int g_node, int NUMNP, int NELEM, double *roe_R, double *roe_a, MNODE *node, int i,\
TRANSFORMED *metric, DETERM *deter, double **f_iplus, double **f_iminus);

