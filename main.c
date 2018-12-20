
/******************SRI ANJINEYA***********************/
/******************OM GAM GAM GANAPATHAYE NAMAHA******/
/***************************OM************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<string.h>
#include<strings.h>
#include"function.h"
#include<mpi.h>
#include<omp.h>

double Epsilon = 10e-8;
double del_zeta=1.0;
double del_eta=1.0;
//0.000071429
double Mach = 3.1;
double Free_t = 273.0;
float Free_rho = 1.225;
float Free_mu = 1.7894e-5;
float Free_pres = 101325.0;
float Free_kine = 1.4;
float Free_sound = 340.28;
float univ_gas = 287.0;
float char_length = 1.0;
//float gamma = 1.4;
double Pr = 0.72;
double del_t = 0.01;
int restart_num = 1;
double Reyl = 220.54874;
int iterations = 20000000;
char *file_ext = "cfx5";             /**********neu if file made in gambit******cfx5 if file made in icemcfd***********/
double back_pressure = 2.74;

int main(int argc, char **argv)
{
	MPI_Init(&argc,&argv );
	int memsize, memsize1;
	char *buffer, *buffer2, *Bcast_buffer1, *Bcast_buffer2, *Bcast_buffer3, *Bcast_buffer4;
	int slaves[1], myrank, size;
	slaves[0]=0;
	MPI_Group group_a,group_b;
	MPI_Comm comm_slaves;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Comm_size(MPI_COMM_WORLD,&size);
  	MPI_Status status;
	//MPI_Request request;
	MPI_Comm_group(MPI_COMM_WORLD,&group_a);
	MPI_Group_excl(group_a,1,slaves,&group_b);
	MPI_Comm_create(MPI_COMM_WORLD,group_b,&comm_slaves);
	
	MPI_Barrier(MPI_COMM_WORLD);

	int i,j, NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL, m, inl, out, wal, bou, NE, IBCODE, all_bou, bound, memory, memory_elem, nodes, no_of_nodes, vn, n[3];
	int out_old, inl_old, wal_old, bou_old, out_mem, out_nl, inl_mem, inl_nl, wal_mem, wal_nl, bou_mem, bou_nl;
	int inl_elem, out_elem, wal_elem, bou_elem;
	int temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	int gar, gar1, gar2, element, elem1, n1, count,k, orient[2], line1, shift, line2, len, len2, shift2, gar3,gar4, gar5, value;
	int P0, P1, P2, P3, P4, P5, P6, P7;
	double angle[4], max_angle, denominator, determinant;
	char garbage1[12];
	int *inlet, *outlet, *wall, *boundary, *boundary_elements, *all_boundary_nodes, *temp_boundary_nodes;
	int *inlet_node, *outlet_node, *wall_node, *boundary_node, g_node, g_elem, all_bou_node;
	int corner_element[20], corner_node[20], cor, wal_node, initial, inl_node, out_node, corn_node, bou_node, temp, restart, position;
	double *det, temp_d;
	double *DUCROS;
//	int icemcfd, gambit;
	double **u, **v, **w, **rho, **p, **t, **mu, *tauzz, *tauxx, *tauzx, *tauex, *tauze, *tauee, *qz, *qe, *qx, **a, **e, **ki;
	//double *u_zeta, *u_eta, *u_xi, *v_zeta, *v_eta, *v_xi, *w_zeta, *w_eta, *w_xi, *t_zeta, *t_eta, *t_xi, *U1_zeta, *U2_zeta, *U3_zeta, *U4_zeta, *U5_zeta, *U1_eta, *U2_eta, *U3_eta, *U4_eta, *U5_eta, *U1_xi, *U2_xi, *U3_xi, *U4_xi, *U5_xi;
	//double eigen_f[5], eigen_E[5], eigen_G[5];
	double Qi_iplus_half[5], Qi_iminus_half[5], Qj_iplus_half[5], Qj_iminus_half[5], Qk_iplus_half[5], Qk_iminus_half[5];
	double  **Qi_iminus, **Qi_iplus, **Qi_iminus_n, **Qi_iplus_n, **Qi_half_p, **Qi_half_np, **Qi_half_m, **Qi_halfn_m; 
	double **Qj_iminus, **Qj_iplus, **Qj_iminus_n, **Qj_iplus_n, **Qj_half_p, **Qj_half_np, **Qj_half_m, **Qj_halfn_m;
	double **Qk_iminus, **Qk_iplus, **Qk_iminus_n, **Qk_iplus_n, **Qk_half_p, **Qk_half_np, **Qk_half_m, **Qk_halfn_m;
	double **IS_Qim, **IS_Qinm, **IS_Qip, **IS_Qinp, **IS_Qjm, **IS_Qjnm,  **IS_Qjp, **IS_Qjnp, **IS_Qkm, **IS_Qknm,  **IS_Qkp, **IS_Qknp, **w_Qip, **w_Qinp, **w_Qim, **w_Qinm, **w_Qjp, **w_Qjnp, **w_Qjm, **w_Qjnm, **w_Qkp, **w_Qknp, **w_Qkm, **w_Qknm, **W_Qip, **W_Qinp, **W_Qim, **W_Qinm, **W_Qjp, **W_Qjnp, **W_Qjm, **W_Qjnm, **W_Qkp, **W_Qknp, **W_Qkm, **W_Qknm;
	double ***U, ***L, **E, **F, **G, **Ev, **Fv, **Gv, **F1, **E1, **G1, **Fv1, **Ev1, **Gv1;
	//double **E_c, **F_c, **G_c;
	double *roe_rho_ip, *roe_u_ip, *roe_v_ip, *roe_w_ip, *roe_h_ip, *roe_rho_im, *roe_u_im, *roe_v_im, *roe_w_im, *roe_h_im, *roe_rho_jp, *roe_u_jp, *roe_v_jp, *roe_w_jp, *roe_h_jp, *roe_rho_jm, *roe_u_jm, *roe_v_jm, *roe_w_jm, *roe_h_jm, *roe_rho_kp, *roe_u_kp, *roe_v_kp, *roe_w_kp, *roe_h_kp, *roe_rho_km, *roe_u_km, *roe_v_km, *roe_w_km, *roe_h_km;
	double **r_eigen_Qip, **r_eigen_Qjp,  **r_eigen_Qkp, **r_eigen_Qim, **r_eigen_Qjm, **r_eigen_Qkm, **l_eigen_Qip, **l_eigen_Qjp, **l_eigen_Qkp, **l_eigen_Qim, **l_eigen_Qjm, **l_eigen_Qkm;
	//double **eigen_val_Qi, **eigen_val_Qj, **eigen_val_Qk, **eigen_i_Qi, **eigen_i_Qj,  **eigen_i_Qk;
	double *dF, *dE, *dG, *dFv, *dEv, *dGv, *final_U;
	double *roe_R, *roe_a, *del_cfl, *v_dash1;
//	double Reyl;
	double *TAU_SGS_XY, *TAU_SGS_YZ, *TAU_SGS_XZ, *TAU_SGS_XX, *TAU_SGS_YY, *TAU_SGS_ZZ, *H_SGS_X, *H_SGS_Y, *H_SGS_Z, *D_SGS_X, *D_SGS_Y, *D_SGS_Z;
	int gambit, icemcfd;
	
	initial = 0;
	MNODE *node;
	singul *singular;
	RESIDUAL **div;
	JACOB *jacobian;
	ELEM *CD; 
	ELEM *wq;
	//MBLOCK *block;
	TRANSFORMED *metric;
	DETERM *deter;
	Qip **Q;
	JAC *jac;
	TEMP_JAC *jaco;
	TMP_NODE *temp_node;
	char line[600];
	char line_n[80];
	char filename[100];
	int itera1, itera2;
	double temp_u, temp_v, temp_w, temp_p, temp_t, temp_e, temp_mu, temp_rho;
	int itera;
	FILE *fp= NULL;
	
	i = 0;
	fp = fopen("anderson_intake_3d.neu","r"); 
	while(fgets(line, 150, fp) != NULL)
	{
		i++;
		if(i == 7)
		{
			sscanf(line, "%d %d %d %d %d %d", &NUMNP, &NELEM, &NGRPS, &NBSETS, &NDFCD, &NDFVL);
			break;
		}		
	}
	fclose(fp);
	itera1 = -1;
	itera2 = -1;
	fp = fopen("restart_file_1.neu","rt");
	if (fp != NULL)
	{
		fgets(line, 1000, fp);
		sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %d", &temp_u, &temp_v, &temp_w, &temp_p, &temp_t, &temp_rho, &temp_e, &temp_mu, &itera);
		//itera++;
	}	
	if (fp != NULL)
	{
		fclose(fp);
	}
	
	itera1 = itera;
	fp = fopen("restart_file_2.neu","rt");
	if (fp != NULL)
	{
		fgets(line, 1000, fp);
		sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %d", &temp_u, &temp_v, &temp_w, &temp_p, &temp_t, &temp_rho, &temp_e, &temp_mu, &itera);
		//itera++;
	}
	if (fp != NULL)
	{
		fclose(fp);
	}
	itera2 = itera;
	
	if (itera1 > itera2 && itera1 >=0 && itera2 >= 0)
	{
		restart_num = 1;
		gar1 = itera1;
	}
	else if (itera2 > itera1 && itera1 >=0 && itera2 >= 0)
	{
		restart_num = 2;
		gar1 = itera2;
	}
	gar = restart_num;
	
	sprintf(filename,"restart_file_%d.neu",restart_num);
	fp = fopen(filename,"rt");
	if (fp != NULL)
	{
		i=1;
		while(fgets(line, 300, fp) != NULL)
		{
			sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %d", &temp_u, &temp_v, &temp_w, &temp_p, &temp_t, &temp_rho, &temp_e, &temp_mu, &itera);
			i++;
			
			if (gar1 != itera)
			{
				i--;
				if (myrank == 0)
				{
					printf("breaking\n");
				}
				break;
			}
			itera = -1;
			 memset( line, 0x00, sizeof(line) );
			line[0] = 0;
			if(i > NUMNP  )
			{
				break;
			} 
		}
		if (i <= NUMNP && gar == 1)
		{
			restart_num++;
		}
		if (i <= NUMNP && gar == 2)
		{
			restart_num--;
		}		
	}	
	if (fp != NULL)
	{
		fclose(fp);
	}
	
	
	if (myrank == 0)
	{		
		printf("restart_num = %d\n", restart_num);
		int restart, iter;
		FILE *fnode;
		FILE *felem;
		    /*****************element file*****************/
		FILE *fp1 = NULL;	/*****************element file*****************/
		FILE *fp2 = NULL;	/*****************node file*****************/
		fp = fopen("anderson_intake_3d.neu","r");    /*****************element file*****************/
		fp1 = fopen("anderson_intake_3d.neu","r");	/*****************element file*****************/
		fp2 = fopen("anderson_intake_3d.neu","r");
		i = 0;

		/******************************************************************************************/
/*
		if (fp == NULL)
		{
			printf("mesh input file not found\n");
			printf("  :(     :(    :(      \n");
			printf("exiting program\n");
			exit(0);
		}
		
		fnode = fopen("node_neighbour.neu","r");
		if (fnode != NULL);
		{
			fclose(fnode);
			system("rm -rf node_neighbour.neu elem_neighbour.neu");
		}
		
		*/
		system("rm -rf node_neighbour.neu elem_neighbour.neu");
		if(strcmp(file_ext,"neu")==0 )
		{
			gambit = 1;
			icemcfd = 0;
		}
		
		if(strcmp(file_ext,"cfx5")==0 )
		{
			gambit = 0;
			icemcfd = 1;
		}
		
		
		/*********************Read Number of node and elements from file***************************/
		//char line[80] = "10";
		//fgets(line, 200, fp);
		while(fgets(line, 200, fp) != NULL)
		{
			i++;
			if(i == 7)
			{
				sscanf(line, "%d %d %d %d %d %d", &NUMNP, &NELEM, &NGRPS, &NBSETS, &NDFCD, &NDFVL);
				printf("Number of nodal points = %d\nNumber of Elements = %d\n", NUMNP, NELEM);
				break;
			}		
		}
		
		/************Obtaining number of boundary nodes*******************************************/
		no_of_nodes = 0;
		inl=0;
		out=0;
		wal=0;
		bou=0;
		out_mem=1;
		out_nl =0;
		out_old=0;
		inl_mem=1;
		inl_nl =0;
		inl_old=0;
		wal_mem=1;
		wal_nl =0;
		wal_old=0;
		bou_mem=1;
		bou_nl =0;
		bou_old=0;
		all_bou=0;
		bound =0;
		cor=0;
		while(fgets(line_n, 200, fp2) != NULL)
		{
			//memset(garbage1,0*sizeof(garbage1));
			//memset(garbage2,0*sizeof(garbage2));
			fscanf(fp2,"%s",garbage1);
			if(strcmp(garbage1,"BOUNDARY")==0)
			{
				//printf("%s\n","found string");
				//printf("%d\n",__LINE__);
				break;
			}
		}
		fgets(line_n, 200, fp2);	

		for(k=0;k<NBSETS;k++)
		{
			while(fgets(line_n, 200, fp2) != NULL)
			{
				//memset(garbage1,0*sizeof(garbage1));
				//memset(garbage2,0*sizeof(garbage2));
				sscanf(line_n,"%s %d %d %d %d",garbage1, &gar1, &NE, &gar3, &gar4);
				//printf("%s\n",line);
			
				if(strncmp(garbage1,"inlet",5)==0)
				{
					no_of_nodes = no_of_nodes+NE;
				}
				
				else if(strncmp(garbage1,"outlet",6)==0)
				{
					no_of_nodes = no_of_nodes+NE;
				}
				
				else if(strncmp(garbage1,"wall",4)==0)
				{
					no_of_nodes = no_of_nodes+NE;
				}
				
				else if(strncmp(garbage1,"boundary",8)==0)
				{
					no_of_nodes = no_of_nodes+NE;
					
				}
				
				if(strncmp(garbage1,"inlet",5)==0)
				{
					inl_old=inl_old+NE;
					IBCODE=1;
					if (inl_mem==inl_nl)
					{
						inlet_node=(int *)realloc(inlet_node,((NE+10)+(inl_old))*sizeof(int));
						
						//printf("%s\n","success");
						inl_mem++;
					}
					else
					{
						//printf("IBCODE %d %d\n",IBCODE,k);
						inlet_node=(int *)malloc((NE+10)*sizeof(int));
						if (inlet_node == NULL)
						{
							printf("LINE %d out of memory\n",__LINE__);
							exit(1);
						}
						
						//printf("%s\n","inlet");
						
					}
					inl_nl++;
				}
				
				else if(strncmp(garbage1,"outlet",6)==0)
				{
					out_old=out_old+NE;
					IBCODE=2;
					if (out_mem==out_nl)
					{
						outlet_node=(int *)realloc(outlet_node,((NE+10)+(out_old))*sizeof(int));
						out_mem++;
					}
					else
					{
						//printf("IBCODE %d %d\n",IBCODE,k);
						outlet_node=(int *)malloc((NE+10)*sizeof(int));
						if (outlet_node == NULL)
						{
							printf("LINE %d out of memory\n",__LINE__);
							exit(1);
						}
					}
					out_nl++;
				}
				
				else if(strncmp(garbage1,"wall",4)==0)
				{
					//wal_mem++;
					wal_old=wal_old+NE;
					IBCODE=3;
					if (wal_mem==wal_nl)
					{
						//printf("%s\n","inside");
						wall_node=(int *)realloc(wall_node,((NE+10)+(wal_old))*sizeof(int));
						wal_mem++;
					}
					else
					{
						//printf("IBCODE %d %d\n",IBCODE,k);
						wall_node=(int *)malloc((NE+10)*sizeof(int));
						if (wall_node == NULL)
						{
							printf("LINE %d out of memory\n",__LINE__);
							exit(1);
						}
					}
					wal_nl++;
				}
				
				else if(strncmp(garbage1,"boundary",8)==0)
				{
				
					bou_old=bou_old+NE;
					IBCODE=4;
					if (bou_mem==bou_nl)
					{
						boundary_node=(int *)realloc(boundary_node,((NE+10)+(bou_old))*sizeof(int));
						bou_mem++;
					}
					else
					{
						boundary_node=(int *)malloc((NE+10)*sizeof(int));
						if (boundary_node == NULL)
						{
							printf("LINE %d out of memory\n",__LINE__);
							exit(1);
						}
					}
					bou_nl++;
				}				
			}
		}

		fclose(fp2);
		nodes = no_of_nodes;

		node = (MNODE*)malloc(((no_of_nodes*6)+(NUMNP+10))*sizeof(MNODE));
		temp_node = (TMP_NODE*)malloc(((no_of_nodes*6)+(NUMNP+10))*sizeof(TMP_NODE));
		CD = (ELEM*)malloc(((no_of_nodes*6)+(NELEM+10))*sizeof(ELEM));
		//wq = (ELEM*)malloc((NELEM+10)*sizeof(ELEM));
		nodes = (no_of_nodes*15)+(NUMNP+1);	
		
		if (node == NULL || temp_node == NULL || CD == NULL )
		{
			printf("LINE %d out of memory\n",__LINE__);
		}
		
/*		memset(node, 0, sizeof(node));
		memset(CD, 0, sizeof(CD));
		memset(wq, 0, sizeof(wq));
*/		
		
		for(i=0; i<=NUMNP; i++)
		{
			node[i].x=0.0;
			node[i].y=0.0;
			node[i].z=0.0;
			
			node[i].e[0]=0;
			node[i].e[1]=0;
			node[i].e[2]=0;
			node[i].e[3]=0;	
			node[i].e[4]=0;
			node[i].e[5]=0;
			node[i].e[6]=0;	
			node[i].e[7]=0;	
			
			temp_node[i].e[0]=0;
			temp_node[i].e[1]=0;
			temp_node[i].e[2]=0;
			temp_node[i].e[3]=0;	
			temp_node[i].e[4]=0;
			temp_node[i].e[5]=0;
			temp_node[i].e[6]=0;	
			temp_node[i].e[7]=0;
			temp_node[i].e[8]=0;
			
			node[i].n_n[0]=0;
			node[i].n_n[1]=0;
			node[i].n_n[2]=0;
			node[i].n_n[3]=0;
			node[i].n_n[4]=0;
			node[i].n_n[5]=0;	
			
			node[i].ID=0;		
		//	node[i].location=0;
		//	node[i].corner=0;
			node[i].corner_ID=0;
			node[i].val=0;
			node[i].loc=0;
		//	node[i].weno=0;
			node[i].proc=0;
			node[i].local=0;
			node[i].global=0;
			node[i].req=0;
		}

		for (i=0; i<= NELEM; i++)
		{
			CD[i].connect[0] = 0;
			CD[i].connect[1] = 0;
			CD[i].connect[2] = 0;
			CD[i].connect[3] = 0;
			CD[i].connect[4] = 0;
			CD[i].connect[5] = 0;
			CD[i].connect[6] = 0;
			CD[i].connect[7] = 0;
		}
		/********************************************************************************************/
		/********************************************************************************************/
		/**********************************************************************************************
				  storing all boundary elements in a single array.... 
	              this array will be freed before solver is executed
		***********************************************************************************************/			  
/*		bound=0;
		cor=0;
		all_bou=inl+out+wal+bou;
		boundary_elements=(int*)malloc((all_bou+10)*sizeof(int));
		memset(boundary_elements,0,(all_bou+10)*sizeof(int));
		for (i=0;i<inl;i++)
		{
			boundary_elements[bound]=inlet[i];
			CD[inlet[i]].type=1;
			bound++;
		}
		for (i=0;i<out;i++)
		{
			boundary_elements[bound]=outlet[i];
			if(CD[outlet[i]].type == 0)
			{
				CD[outlet[i]].type=2;
			}
			else
			{
				corner_element[cor]=outlet[i];
				cor++;
			}
			bound++;
		}
		
		for (i=0;i<wal;i++)
		{
			boundary_elements[bound]=wall[i];
			if (CD[wall[i]].type == 0)
			{
				CD[wall[i]].type=3;
			}
			else
			{
				corner_element[cor]=wall[i];
				cor++;
			}
			bound++;
		}
		for (i=0;i<bou;i++)
		{
			boundary_elements[bound]=boundary[i];
			if (CD[boundary[i]].type == 0)
			{
				CD[boundary[i]].type=4;
			}
			else
			{
				corner_element[cor]=boundary[i];
				cor++;
			}
			bound++;
		}
*/
		/***********************Storing x and y coordinates of each node*****************************/
		i = 1;
		j = 0;
		while(fgets(line, 200, fp) != NULL)
		{
			if(j > 1 && j<=NUMNP+2)
			{
				sscanf(line, "%d %lf %lf %lf", &gar, &node[i].x, &node[i].y, &node[i].z);
				node[i].ID=0;
				i++;	
			}
			j++;
			if (j>NUMNP+1)
			{
				break;
			}
		}
		
		/*********************************************************************************************/

		/*****************************storing node numbers of each element****************************/	
		i = 1;
		j = 0;
		while(fgets(line, 200, fp) != NULL)
		{
			if (j>1 && j<=NELEM+1)
			{
				if (icemcfd == 1)
				{
					sscanf(line, "%d %d %d %d %d %d %d %d %d %d", &gar, &gar1, &gar2, &CD[i].connect[4], &CD[i].connect[5], &CD[i].connect[7], &CD[i].connect[6], &CD[i].connect[0], &CD[i].connect[1], &CD[i].connect[3]);
					fgets(line, 150, fp);
					sscanf(line, "%d", &CD[i].connect[2]);
				}
				if (gambit == 1)
				{
					sscanf(line, "%d %d %d %d %d %d %d %d %d %d", &gar, &gar1, &gar2, &CD[i].connect[0], &CD[i].connect[1], &CD[i].connect[3], &CD[i].connect[2], &CD[i].connect[4], &CD[i].connect[5], &CD[i].connect[7]);
					fgets(line, 150, fp);
					sscanf(line, "%d", &CD[i].connect[6]);
				}
				
				temp_node[CD[i].connect[0]].e[node[CD[i].connect[0]].val++] = i;
				temp_node[CD[i].connect[1]].e[node[CD[i].connect[1]].val++] = i;
				temp_node[CD[i].connect[2]].e[node[CD[i].connect[2]].val++] = i;
				temp_node[CD[i].connect[3]].e[node[CD[i].connect[3]].val++] = i;
				temp_node[CD[i].connect[4]].e[node[CD[i].connect[4]].val++] = i;
				temp_node[CD[i].connect[5]].e[node[CD[i].connect[5]].val++] = i;
				temp_node[CD[i].connect[6]].e[node[CD[i].connect[6]].val++] = i;
				temp_node[CD[i].connect[7]].e[node[CD[i].connect[7]].val++] = i;
				
				i++;
			}
			if (j == NELEM+1)
			{
				break;
			}
			j++;
		}
		
		fclose(fp);
		/**********************************************************************************************/
		/******************************************Reading boundary nodes from file*************************************/
		fp2 = fopen("anderson_intake_3d.neu","r");
		inl=0;
		out=0;
		wal=0;
		bou=0;
		out_mem=1;
		out_nl =0;
		out_old=0;
		inl_mem=1;
		inl_nl =0;
		inl_old=0;
		wal_mem=1;
		wal_nl =0;
		wal_old=0;
		bou_mem=1;
		bou_nl =0;
		bou_old=0;
		all_bou=0;
		bound =0;
		cor=0;
		while(fgets(line_n, 200, fp2) != NULL)
		{
			//memset(garbage1,0*sizeof(garbage1));
			//memset(garbage2,0*sizeof(garbage2));
			fscanf(fp2,"%s",garbage1);
			if(strcmp(garbage1,"BOUNDARY")==0)
			{
				//printf("%s\n","found string");
				//printf("%d\n",__LINE__);
				break;
			}
		}
		fgets(line_n, 200, fp2);	

		for(k=0;k<NBSETS;k++)
		{
			while(fgets(line_n, 200, fp2) != NULL)
			{
				sscanf(line_n,"%s %d %d %d %d",garbage1, &gar1, &NE, &gar3, &gar4);
			
				if(strncmp(garbage1,"inlet",5)==0)
				{
					inl_old=inl_old+NE;
					IBCODE=1;
				}
				
				else if(strncmp(garbage1,"outlet",6)==0)
				{
					out_old=out_old+NE;
					IBCODE=2;
				}
				
				else if(strncmp(garbage1,"wall",4)==0)
				{
					wal_old=wal_old+NE;
					IBCODE=3;
				}
				
				else if(strncmp(garbage1,"boundary",8)==0)
				{
					bou_old=bou_old+NE;
					IBCODE=4;
				}
				
				switch (IBCODE)
				{
					case 1:
						gar5=0;
						while(fgets(line_n, 200, fp2) != NULL)
						{
							if(gar5==NE)
							{
								//printf("%s\n","inlet EOS");
								break;
							}
							sscanf(line_n,"%d ",&inlet_node[inl]);
							if (node[inlet_node[inl]].ID == 0 )
							{
								node[inlet_node[inl]].ID=1;
							}
							else
							{
								if (node[inlet_node[inl]].ID != 3 )
								{
									node[inlet_node[inl]].ID=1;
								}
								//corner_node[cor]=inlet_node[inl];
//								node[inlet_node[inl]].corner=1;
							}
							inl++;
							gar5++;	 
						}
						fgets(line_n,  200, fp2);
					break;
			
					case 2:
						gar5=0;
						while(fgets(line_n, 200, fp2) != NULL)
						{
							if (gar5==NE)
							{
								break;
							}
							sscanf(line_n,"%d ",&outlet_node[out]);
							if (node[outlet_node[out]].ID == 0)
							{
								node[outlet_node[out]].ID=2;
							}
							else
							{
								//corner_node[cor]=outlet_node[out];
//								node[outlet_node[out]].corner=1;
								cor++;
							}
							out++;					
							gar5++;	
						}
						fgets(line_n, 200, fp2);		
					break;
			
					case 3:
						gar5=0;
						while(fgets(line_n, 200, fp2) != NULL)
						{
							if(gar5==NE)
							{
								break;
							}
							sscanf(line_n,"%d ",&wall_node[wal]);
							if (node[wall_node[wal]].ID == 0)
							{
								node[wall_node[wal]].ID=3;
							}
							else
							{
								node[wall_node[wal]].ID=3;
								//corner_node[cor]=wall_node[wal];
//								node[wall_node[wal]].corner=1;
								cor++;
							}
							wal++;					
							gar5++;	
						}
						fgets(line_n, 200, fp2);		
					break;
			
					case 4:
						gar5=0;
						while(fgets(line_n, 200, fp2) != NULL)
						{
							if(gar5==NE)
							{
								break;
							}
							boundary_node[bou] = 0;
							sscanf(line_n,"%d",&boundary_node[bou]);
							if (node[boundary_node[bou]].ID == 0)
							{
								node[boundary_node[bou]].ID=4;
							}
							else
							{
								//corner_node[cor]=boundary_node[bou];
//								node[boundary_node[bou]].corner=1;
								cor++;
							}
							bou++;					
							gar5++;	
						}
						fgets(line_n, 200, fp2);		
					break;
				}
				IBCODE =10;
				break;
			}
		}
		fclose(fp2);
		corn_node = cor;
		wal_node = wal;	
		inl_node = inl;
		out_node = out;
		bou_node = bou;

		/*************************Looking for node_NEIGHBOUR.neu FILE******************************************/
		printf("started loking for node_neighbour.neu\n");
		fnode =fopen("node_neighbour.neu","rt"); 
		if(fnode == NULL)
		{
			/**************************Rearranging node numbers of each element****************************/
			printf("****************************************************\n");
			printf("node_neighbour.neu not found\n");
			printf("creating node_neigh file\n");
			printf("****************************************************\n");
//			#pragma omp parallel private(max_angle, line1, line2, len, len2, denominator, angle, shift, shift2, i) shared(node, CD)
//			{	
//				#pragma omp for
		//		if (gambit == 1)
				{	
					for(element = 1; element <=NELEM; element++)
					{
						temp = 0;
						if (fabs(node[CD[element].connect[0]].x-node[CD[element].connect[4]].x) >= fabs(node[CD[element].connect[0]].y-node[CD[element].connect[4]].y) && fabs(node[CD[element].connect[0]].x-node[CD[element].connect[4]].x) >= fabs(node[CD[element].connect[0]].z-node[CD[element].connect[4]].z))
						{
							value = 1;
						}
						if (fabs(node[CD[element].connect[0]].y-node[CD[element].connect[4]].y) >= fabs(node[CD[element].connect[0]].x-node[CD[element].connect[4]].x) && fabs(node[CD[element].connect[0]].y-node[CD[element].connect[4]].y) >= fabs(node[CD[element].connect[0]].z-node[CD[element].connect[4]].z))
						{
							value = 2;
						}
						if (fabs(node[CD[element].connect[0]].z-node[CD[element].connect[4]].z) >= fabs(node[CD[element].connect[0]].x-node[CD[element].connect[4]].x) && fabs(node[CD[element].connect[0]].z-node[CD[element].connect[4]].z) >= fabs(node[CD[element].connect[0]].y-node[CD[element].connect[4]].y))
						{
							value = 3;
						}
				
						switch (value)
						{
							case 1:
							
								max_angle = 0;
								line1 = 0;
								line2 = 0;
								len = 0;
								len2= 0;
								for (i = 0; i<4; i++)
								{
									if (i<3)
									{
										denominator =(node[CD[element].connect[i+1]].z-node[CD[element].connect[i]].z);
										if (denominator <1e-16 && denominator >= 0.0)
										{
											angle[i]=90;
										}
										else
										{
											angle[i]=atan(fabs(((node[CD[element].connect[i+1]].y-node[CD[element].connect[i]].y)/ denominator)));
										}
									}
									else
									{
										denominator = (node[CD[element].connect[i]].z-node[CD[element].connect[0]].z);
										if (denominator <1e-16 && denominator >= 0.0)
										{
											angle[i]=90;
										}
										else
										{
											angle[i]=atan(fabs(((node[CD[element].connect[i]].y-node[CD[element].connect[0]].y)/denominator)));
										}
									}
									if (angle[i] >max_angle)
									{
										max_angle = angle[i];
										line1 = i;
									}
								}
								if (line1 == 0 ) 
								{
									line2 = 2;
								}
								else if (line1 == 1)
								{
									line2 = 3;
								}
								else if (line1 == 2)
								{
									line2 = 0;
								}
								else
								{
									line2 = 1;
								}
								if (line1 == 0 )
								{
									if (node[CD[element].connect[line1]].y > node[CD[element].connect[line1+1]].y) 
									{
										if (node[CD[element].connect[0]].z > node[CD[element].connect[line2+1]].z)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											shift2 = CD[element].connect[2];
											CD[element].connect[2] = CD[element].connect[1];
											CD[element].connect[1] = shift;
											CD[element].connect[3] = shift2;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											shift2 = CD[element].connect[6];
											CD[element].connect[6] = CD[element].connect[5];
											CD[element].connect[5] = shift;
											CD[element].connect[7] = shift2;					
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[1];
											CD[element].connect[1] = shift;
										
											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[5];
											CD[element].connect[5] = shift;
										}
									}
									else if (node[CD[element].connect[0]].z > node[CD[element].connect[line2+1]].z)
									{
										shift = CD[element].connect[2];
										CD[element].connect[2] = CD[element].connect[0];
										CD[element].connect[0] = shift;
									
										shift = CD[element].connect[6];
										CD[element].connect[6] = CD[element].connect[4];
										CD[element].connect[4] = shift;
									}
									else
									{
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										shift2 = CD[element].connect[3];
										CD[element].connect[3] = shift;
										CD[element].connect[1] = CD[element].connect[2];
										CD[element].connect[2] = shift2;	

										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										shift2 = CD[element].connect[7];
										CD[element].connect[7] = shift;
										CD[element].connect[5] = CD[element].connect[6];
										CD[element].connect[6] = shift2;
									}
								}
				
								if (line1 == 2)
								{
									if (node[CD[element].connect[line1]].y < node[CD[element].connect[line1+1]].y) 
									{
										if (node[CD[element].connect[line1]].z < node[CD[element].connect[1]].z)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											shift2 = CD[element].connect[2];
											CD[element].connect[2] = CD[element].connect[1];
											CD[element].connect[1] = shift;
											CD[element].connect[3] = shift2;	

											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											shift2 = CD[element].connect[6];
											CD[element].connect[6] = CD[element].connect[5];
											CD[element].connect[5] = shift;
											CD[element].connect[7] = shift2;
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[1];
											CD[element].connect[1] = shift;		

											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[5];
											CD[element].connect[5] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].z < node[CD[element].connect[1]].z)
									{
										shift = CD[element].connect[2];
										CD[element].connect[2] = CD[element].connect[0];
										CD[element].connect[0] = shift;	

										shift = CD[element].connect[6];
										CD[element].connect[6] = CD[element].connect[4];
										CD[element].connect[4] = shift;
									}
									else
									{
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										shift2 = CD[element].connect[3];
										CD[element].connect[3] = shift;
										CD[element].connect[1] = CD[element].connect[2];
										CD[element].connect[2] = shift2;		

										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										shift2 = CD[element].connect[7];
										CD[element].connect[7] = shift;
										CD[element].connect[5] = CD[element].connect[6];
										CD[element].connect[6] = shift2;
									}
								}
				
								if (line1 == 1) 
								{
									len = 2;
									len2 = 0;
						
									if (node[CD[element].connect[line1]].y < node[CD[element].connect[len]].y) 
									{
										if (node[CD[element].connect[line1]].z > node[CD[element].connect[len2]].z)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											CD[element].connect[3] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[2];
											CD[element].connect[2] = shift;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											CD[element].connect[7] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[6];
											CD[element].connect[6] = shift;
										}
										else
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[2];
											CD[element].connect[2] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[3];
											CD[element].connect[3] = shift;		

											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[6];
											CD[element].connect[6] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[7];
											CD[element].connect[7] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].z < node[CD[element].connect[len2]].z)
									{
										shift = CD[element].connect[3];
										CD[element].connect[3] = CD[element].connect[2];
										CD[element].connect[2] = shift;
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										CD[element].connect[1] = shift;	

										shift = CD[element].connect[7];
										CD[element].connect[7] = CD[element].connect[6];
										CD[element].connect[6] = shift;
										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										CD[element].connect[5] = shift;
									}
								}
								if (line1 == 3)
								{
									if (line1 == 3)
									{
										len = 0;
										len2 = 2;
									}

									if (node[CD[element].connect[line1]].y > node[CD[element].connect[len]].y) 
									{
										if (node[CD[element].connect[line1]].z > node[CD[element].connect[len2]].z)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[2];
											CD[element].connect[2] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[3];
											CD[element].connect[3] = shift;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[6];
											CD[element].connect[6] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[7];
											CD[element].connect[7] = shift;
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[0];
											CD[element].connect[0] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[2];
											CD[element].connect[2] = shift;		
										
											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[4];
											CD[element].connect[4] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[6];
											CD[element].connect[6] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].z > node[CD[element].connect[len2]].z)
									{
										shift = CD[element].connect[3];
										CD[element].connect[3] = CD[element].connect[2];
										CD[element].connect[2] = shift;
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										CD[element].connect[1] = shift;	

										shift = CD[element].connect[7];
										CD[element].connect[7] = CD[element].connect[6];
										CD[element].connect[6] = shift;
										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										CD[element].connect[5] = shift;
									}
								}					
								P0 = CD[element].connect[0];
								P1 = CD[element].connect[1];
								P2 = CD[element].connect[2];
								P3 = CD[element].connect[3];
								P4 = CD[element].connect[4];
								P5 = CD[element].connect[5];
								P6 = CD[element].connect[6];
								P7 = CD[element].connect[7];
							
								if (node[CD[element].connect[0]].x-node[CD[element].connect[4]].x <= 0.0 && temp == 0)
								{
									CD[element].connect[0] = P2;
									CD[element].connect[1] = P6;
									CD[element].connect[2] = P5;
									CD[element].connect[3] = P1;
									CD[element].connect[4] = P3;
									CD[element].connect[5] = P7;
									CD[element].connect[6] = P4;
									CD[element].connect[7] = P0;	
								
									temp = 1;
								}
								if (node[CD[element].connect[0]].x-node[CD[element].connect[4]].x >= 0.0 && temp == 0)
								{
									CD[element].connect[0] = P6;
									CD[element].connect[1] = P2;
									CD[element].connect[2] = P1;
									CD[element].connect[3] = P5;
									CD[element].connect[4] = P7;
									CD[element].connect[5] = P3;
									CD[element].connect[6] = P0;
									CD[element].connect[7] = P4;	

									temp = 1;
								}
							break;
						
							case 2:
								max_angle = 0;
								line1 = 0;
								line2 = 0;
								len = 0;
								len2= 0;
								for (i = 0; i<4; i++)
								{
									if (i<3)
									{
										denominator =(node[CD[element].connect[i+1]].x-node[CD[element].connect[i]].x);
										if (denominator <1e-16 && denominator >= 0.0)
										{
											angle[i]=90;
										}
										else
										{
											angle[i]=atan(fabs(((node[CD[element].connect[i+1]].z-node[CD[element].connect[i]].z)/ denominator)));
										}
									}
									else
									{
										denominator = (node[CD[element].connect[i]].x-node[CD[element].connect[0]].x);
										if (denominator <1e-16 && denominator >= 0.0)
										{
											angle[i]=90;
										}
										else
										{
											angle[i]=atan(fabs(((node[CD[element].connect[i]].z-node[CD[element].connect[0]].z)/denominator)));
										}
									}
									if (angle[i] >max_angle)
									{
										max_angle = angle[i];
										line1 = i;
									}
								}
								if (line1 == 0 ) 
								{
									line2 = 2;
								}
								else if (line1 == 1)
								{
									line2 = 3;
								}
								else if (line1 == 2)
								{
									line2 = 0;
								}
								else
								{
									line2 = 1;
								}
								if (line1 == 0 )
								{
									if (node[CD[element].connect[line1]].z > node[CD[element].connect[line1+1]].z) 
									{
										if (node[CD[element].connect[0]].x > node[CD[element].connect[line2+1]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											shift2 = CD[element].connect[2];
											CD[element].connect[2] = CD[element].connect[1];
											CD[element].connect[1] = shift;
											CD[element].connect[3] = shift2;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											shift2 = CD[element].connect[6];
											CD[element].connect[6] = CD[element].connect[5];
											CD[element].connect[5] = shift;
											CD[element].connect[7] = shift2;					
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[1];
											CD[element].connect[1] = shift;
										
											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[5];
											CD[element].connect[5] = shift;
										}
									}
									else if (node[CD[element].connect[0]].x > node[CD[element].connect[line2+1]].x)
									{
										shift = CD[element].connect[2];
										CD[element].connect[2] = CD[element].connect[0];
										CD[element].connect[0] = shift;
									
										shift = CD[element].connect[6];
										CD[element].connect[6] = CD[element].connect[4];
										CD[element].connect[4] = shift;
									}
									else
									{
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										shift2 = CD[element].connect[3];
										CD[element].connect[3] = shift;
										CD[element].connect[1] = CD[element].connect[2];
										CD[element].connect[2] = shift2;	

										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										shift2 = CD[element].connect[7];
										CD[element].connect[7] = shift;
										CD[element].connect[5] = CD[element].connect[6];
										CD[element].connect[6] = shift2;
									}
								}
				
								if (line1 == 2)
								{
									if (node[CD[element].connect[line1]].z < node[CD[element].connect[line1+1]].z) 
									{
										if (node[CD[element].connect[line1]].x < node[CD[element].connect[1]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											shift2 = CD[element].connect[2];
											CD[element].connect[2] = CD[element].connect[1];
											CD[element].connect[1] = shift;
											CD[element].connect[3] = shift2;	

											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											shift2 = CD[element].connect[6];
											CD[element].connect[6] = CD[element].connect[5];
											CD[element].connect[5] = shift;
											CD[element].connect[7] = shift2;
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[1];
											CD[element].connect[1] = shift;		

											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[5];
											CD[element].connect[5] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].x < node[CD[element].connect[1]].x)
									{
										shift = CD[element].connect[2];
										CD[element].connect[2] = CD[element].connect[0];
										CD[element].connect[0] = shift;	

										shift = CD[element].connect[6];
										CD[element].connect[6] = CD[element].connect[4];
										CD[element].connect[4] = shift;
									}
									else
									{
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										shift2 = CD[element].connect[3];
										CD[element].connect[3] = shift;
										CD[element].connect[1] = CD[element].connect[2];
										CD[element].connect[2] = shift2;		

										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										shift2 = CD[element].connect[7];
										CD[element].connect[7] = shift;
										CD[element].connect[5] = CD[element].connect[6];
										CD[element].connect[6] = shift2;
									}
								}
				
								if (line1 == 1) 
								{
									len = 2;
									len2 = 0;
						
									if (node[CD[element].connect[line1]].z < node[CD[element].connect[len]].z) 
									{
										if (node[CD[element].connect[line1]].x > node[CD[element].connect[len2]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											CD[element].connect[3] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[2];
											CD[element].connect[2] = shift;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											CD[element].connect[7] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[6];
											CD[element].connect[6] = shift;
										}
										else
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[2];
											CD[element].connect[2] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[3];
											CD[element].connect[3] = shift;		

											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[6];
											CD[element].connect[6] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[7];
											CD[element].connect[7] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].x < node[CD[element].connect[len2]].x)
									{
										shift = CD[element].connect[3];
										CD[element].connect[3] = CD[element].connect[2];
										CD[element].connect[2] = shift;
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										CD[element].connect[1] = shift;	

										shift = CD[element].connect[7];
										CD[element].connect[7] = CD[element].connect[6];
										CD[element].connect[6] = shift;
										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										CD[element].connect[5] = shift;
									}
								}
								if (line1 == 3)
								{
									if (line1 == 3)
									{
										len = 0;
										len2 = 2;
									}

									if (node[CD[element].connect[line1]].z > node[CD[element].connect[len]].z) 
									{
										if (node[CD[element].connect[line1]].x > node[CD[element].connect[len2]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[2];
											CD[element].connect[2] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[3];
											CD[element].connect[3] = shift;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[6];
											CD[element].connect[6] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[7];
											CD[element].connect[7] = shift;
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[0];
											CD[element].connect[0] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[2];
											CD[element].connect[2] = shift;		
										
											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[4];
											CD[element].connect[4] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[6];
											CD[element].connect[6] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].x > node[CD[element].connect[len2]].x)
									{
										shift = CD[element].connect[3];
										CD[element].connect[3] = CD[element].connect[2];
										CD[element].connect[2] = shift;
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										CD[element].connect[1] = shift;	

										shift = CD[element].connect[7];
										CD[element].connect[7] = CD[element].connect[6];
										CD[element].connect[6] = shift;
										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										CD[element].connect[5] = shift;
									}
								}
													
								P0 = CD[element].connect[0];
								P1 = CD[element].connect[1];
								P2 = CD[element].connect[2];
								P3 = CD[element].connect[3];
								P4 = CD[element].connect[4];
								P5 = CD[element].connect[5];
								P6 = CD[element].connect[6];
								P7 = CD[element].connect[7];
							
								if (node[CD[element].connect[0]].y-node[CD[element].connect[4]].y <= 0.0 && temp == 0)
								{
									CD[element].connect[0] = P0;
									CD[element].connect[1] = P1;
									CD[element].connect[2] = P5;
									CD[element].connect[3] = P4;
									CD[element].connect[4] = P3;
									CD[element].connect[5] = P2;
									CD[element].connect[6] = P6;
									CD[element].connect[7] = P7;    
									temp = 1;
								}
								if (node[CD[element].connect[0]].y-node[CD[element].connect[4]].y >= 0.0 && temp == 0)
								{
									CD[element].connect[0] = P4;
									CD[element].connect[1] = P5;
									CD[element].connect[2] = P1;
									CD[element].connect[3] = P0;
									CD[element].connect[4] = P7;
									CD[element].connect[5] = P6;
									CD[element].connect[6] = P2;
									CD[element].connect[7] = P3; 
									temp = 1;
								}	
							break;					
					
							case 3:
								max_angle = 0;
								line1 = 0;
								line2 = 0;
								len = 0;
								len2= 0;
								for (i = 0; i<4; i++)
								{
									if (i<3)
									{
										denominator =(node[CD[element].connect[i+1]].x-node[CD[element].connect[i]].x);
										if (denominator <1e-16 && denominator >= 0.0)
										{
											angle[i]=90;
										}
										else
										{
											angle[i]=atan(fabs(((node[CD[element].connect[i+1]].y-node[CD[element].connect[i]].y)/ denominator)));
										}
									}
									else
									{
										denominator = (node[CD[element].connect[i]].x-node[CD[element].connect[0]].x);
										if (denominator <1e-16 && denominator >= 0.0)
										{
											angle[i]=90;
										}
										else
										{
											angle[i]=atan(fabs(((node[CD[element].connect[i]].y-node[CD[element].connect[0]].y)/denominator)));
										}
									}
									if (angle[i] >max_angle)
									{
										max_angle = angle[i];
										line1 = i;
									}
								}
								if (line1 == 0 ) 
								{
									line2 = 2;
								}
								else if (line1 == 1)
								{
									line2 = 3;
								}
								else if (line1 == 2)
								{
									line2 = 0;
								}
								else
								{
									line2 = 1;
								}
								if (line1 == 0 )
								{
									if (node[CD[element].connect[line1]].y > node[CD[element].connect[line1+1]].y) 
									{
										if (node[CD[element].connect[0]].x > node[CD[element].connect[line2+1]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											shift2 = CD[element].connect[2];
											CD[element].connect[2] = CD[element].connect[1];
											CD[element].connect[1] = shift;
											CD[element].connect[3] = shift2;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											shift2 = CD[element].connect[6];
											CD[element].connect[6] = CD[element].connect[5];
											CD[element].connect[5] = shift;
											CD[element].connect[7] = shift2;					
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[1];
											CD[element].connect[1] = shift;
										
											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[5];
											CD[element].connect[5] = shift;
										}
									}
									else if (node[CD[element].connect[0]].x > node[CD[element].connect[line2+1]].x)
									{
										shift = CD[element].connect[2];
										CD[element].connect[2] = CD[element].connect[0];
										CD[element].connect[0] = shift;
									
										shift = CD[element].connect[6];
										CD[element].connect[6] = CD[element].connect[4];
										CD[element].connect[4] = shift;
									}
									else
									{
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										shift2 = CD[element].connect[3];
										CD[element].connect[3] = shift;
										CD[element].connect[1] = CD[element].connect[2];
										CD[element].connect[2] = shift2;	

										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										shift2 = CD[element].connect[7];
										CD[element].connect[7] = shift;
										CD[element].connect[5] = CD[element].connect[6];
										CD[element].connect[6] = shift2;
									}
								}
				
								if (line1 == 2)
								{
									if (node[CD[element].connect[line1]].y < node[CD[element].connect[line1+1]].y) 
									{
										if (node[CD[element].connect[line1]].x < node[CD[element].connect[1]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											shift2 = CD[element].connect[2];
											CD[element].connect[2] = CD[element].connect[1];
											CD[element].connect[1] = shift;
											CD[element].connect[3] = shift2;	

											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											shift2 = CD[element].connect[6];
											CD[element].connect[6] = CD[element].connect[5];
											CD[element].connect[5] = shift;
											CD[element].connect[7] = shift2;
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[1];
											CD[element].connect[1] = shift;		

											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[5];
											CD[element].connect[5] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].x < node[CD[element].connect[1]].x)
									{
										shift = CD[element].connect[2];
										CD[element].connect[2] = CD[element].connect[0];
										CD[element].connect[0] = shift;	

										shift = CD[element].connect[6];
										CD[element].connect[6] = CD[element].connect[4];
										CD[element].connect[4] = shift;
									}
									else
									{
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										shift2 = CD[element].connect[3];
										CD[element].connect[3] = shift;
										CD[element].connect[1] = CD[element].connect[2];
										CD[element].connect[2] = shift2;		

										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										shift2 = CD[element].connect[7];
										CD[element].connect[7] = shift;
										CD[element].connect[5] = CD[element].connect[6];
										CD[element].connect[6] = shift2;
									}
								}
				
								if (line1 == 1) 
								{
									len = 2;
									len2 = 0;
						
									if (node[CD[element].connect[line1]].y < node[CD[element].connect[len]].y) 
									{
										if (node[CD[element].connect[line1]].x > node[CD[element].connect[len2]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[3];
											CD[element].connect[3] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[2];
											CD[element].connect[2] = shift;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[7];
											CD[element].connect[7] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[6];
											CD[element].connect[6] = shift;
										}
										else
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[2];
											CD[element].connect[2] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[3];
											CD[element].connect[3] = shift;		

											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[6];
											CD[element].connect[6] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[7];
											CD[element].connect[7] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].x < node[CD[element].connect[len2]].x)
									{
										shift = CD[element].connect[3];
										CD[element].connect[3] = CD[element].connect[2];
										CD[element].connect[2] = shift;
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										CD[element].connect[1] = shift;	

										shift = CD[element].connect[7];
										CD[element].connect[7] = CD[element].connect[6];
										CD[element].connect[6] = shift;
										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										CD[element].connect[5] = shift;
									}
								}
								if (line1 == 3)
								{
									if (line1 == 3)
									{
										len = 0;
										len2 = 2;
									}

									if (node[CD[element].connect[line1]].y > node[CD[element].connect[len]].y) 
									{
										if (node[CD[element].connect[line1]].x > node[CD[element].connect[len2]].x)
										{
											shift = CD[element].connect[0];
											CD[element].connect[0] = CD[element].connect[2];
											CD[element].connect[2] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[3];
											CD[element].connect[3] = shift;
										
											shift = CD[element].connect[4];
											CD[element].connect[4] = CD[element].connect[6];
											CD[element].connect[6] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[7];
											CD[element].connect[7] = shift;
										}
										else
										{
											shift = CD[element].connect[3];
											CD[element].connect[3] = CD[element].connect[0];
											CD[element].connect[0] = shift;
											shift = CD[element].connect[1];
											CD[element].connect[1] = CD[element].connect[2];
											CD[element].connect[2] = shift;		
										
											shift = CD[element].connect[7];
											CD[element].connect[7] = CD[element].connect[4];
											CD[element].connect[4] = shift;
											shift = CD[element].connect[5];
											CD[element].connect[5] = CD[element].connect[6];
											CD[element].connect[6] = shift;	
										}
									}
									else if (node[CD[element].connect[line1]].x > node[CD[element].connect[len2]].x)
									{
										shift = CD[element].connect[3];
										CD[element].connect[3] = CD[element].connect[2];
										CD[element].connect[2] = shift;
										shift = CD[element].connect[0];
										CD[element].connect[0] = CD[element].connect[1];
										CD[element].connect[1] = shift;	

										shift = CD[element].connect[7];
										CD[element].connect[7] = CD[element].connect[6];
										CD[element].connect[6] = shift;
										shift = CD[element].connect[4];
										CD[element].connect[4] = CD[element].connect[5];
										CD[element].connect[5] = shift;
									}
								}
							
								P0 = CD[element].connect[0];
								P1 = CD[element].connect[1];
								P2 = CD[element].connect[2];
								P3 = CD[element].connect[3];
								P4 = CD[element].connect[4];
								P5 = CD[element].connect[5];
								P6 = CD[element].connect[6];
								P7 = CD[element].connect[7];
							
								if (node[CD[element].connect[0]].z-node[CD[element].connect[4]].z <= 0.0 && temp == 0)
								{
									CD[element].connect[0] = P7;
									CD[element].connect[1] = P6;
									CD[element].connect[2] = P5;
									CD[element].connect[3] = P4;
									CD[element].connect[4] = P3;
									CD[element].connect[5] = P2;
									CD[element].connect[6] = P1;
									CD[element].connect[7] = P0;     
									temp = 1;
								}
								if (node[CD[element].connect[0]].z-node[CD[element].connect[4]].z >= 0.0 && temp == 0)
								{
									CD[element].connect[0] = P3;
									CD[element].connect[1] = P2;
									CD[element].connect[2] = P1;
									CD[element].connect[3] = P0;
									CD[element].connect[4] = P7;
									CD[element].connect[5] = P6;
									CD[element].connect[6] = P5;
									CD[element].connect[7] = P4; 
									temp = 1;
								}
							break;				
						}
					//printf("%d -- %d %d %d %d %d %d %d %d\n",element, CD[element].connect[0], CD[element].connect[1], CD[element].connect[2], CD[element].connect[3], CD[element].connect[4], CD[element].connect[5], CD[element].connect[6], CD[element].connect[7]);
					}
				}
			//}	
			/***********************************************************************************************/
			/*****************************Finding Elements connected to each node***************************/
//			#pragma omp parallel private(j,element, n1) shared(node, CD)
//			{	
//				#pragma omp for
				for(n1 = 1; n1 <= NUMNP; n1++)
				{
					for(element = 0; element <= 7; element++)
					{			
						for(j = 0; j < 8; j++)
						{
							if(n1 == CD[temp_node[n1].e[element]].connect[j] )
							{
								if (j == 0)
								{
									node[n1].e[2] = temp_node[n1].e[element];	
								}
								else if (j==1)
								{
									node[n1].e[1] = temp_node[n1].e[element];	
								}
								else if (j==2)
								{
									node[n1].e[5] = temp_node[n1].e[element];	
								}
								else if (j==3)
								{
									node[n1].e[6] = temp_node[n1].e[element];	
								}
								else if (j==4)
								{
									node[n1].e[3] = temp_node[n1].e[element];	
								}
								else if (j==5)
								{
									node[n1].e[0] = temp_node[n1].e[element];	
								}
								else if (j==6)
								{
									node[n1].e[4] = temp_node[n1].e[element];	
								}
								else if (j==7)
								{
									node[n1].e[7] = temp_node[n1].e[element];	
								}
							}
						}
					}
				}
//			}
		

			free (temp_node);
			temp_node = NULL;
			/***********************************************************************************************/
			/*****************	Finding location of node ***************************************************/
			corn_node = 0;
			for (i=1; i<=NUMNP; i++)
			{
				node[i].loc = 0;
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 1;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 2;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 3;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 4;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 5;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 6;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 7;
					node[i].val = node[i].loc;
					node[i].corner_ID = 1;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 8;
					node[i].val = node[i].loc;
					node[i].corner_ID = 2;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 9;
					node[i].val = node[i].loc;
					node[i].corner_ID = 3;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 10;
					node[i].val = node[i].loc;
					node[i].corner_ID = 4;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 11;
					node[i].val = node[i].loc;
					node[i].corner_ID = 5;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 12;
					node[i].val = node[i].loc;
					node[i].corner_ID = 6;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 13;
					node[i].val = node[i].loc;
					node[i].corner_ID = 7;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 14;
					node[i].val = node[i].loc;
					node[i].corner_ID = 8;
					corner_node[corn_node] = i;
					corn_node++;
				}
				if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 15;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 16;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 17;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 18;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 19;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 20;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 21;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 22;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 23;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 24;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 25;
					node[i].val = node[i].loc;
				}
				if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 26;
					node[i].val = node[i].loc;
				}
				
			/*	if((node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0 && (CD[node[i].e[1]].connect[5] != CD[node[i].e[2]].connect[4]) && (CD[node[i].e[1]].connect[6] != CD[node[i].e[2]].connect[7])))
				{
					node[i].loc = 27;
				}
				if((node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0 && (CD[node[i].e[2]].connect[0] != CD[node[i].e[6]].connect[3]) && (CD[node[i].e[3]].connect[4] != CD[node[i].e[7]].connect[7])))
				{
					node[i].loc = 28;
				} */
			}		
			
			for (i = 1; i<= NUMNP; i++)
			{
				/*if(node[CD[node[i].e[5]].connect[1]].loc == 27 && node[CD[node[i].e[6]].connect[1]].loc == 28)
				{
					node[i].loc = 29;
					node[i].ID = 10;
				}
				if(node[CD[node[i].e[5]].connect[3]].loc == 28 && node[CD[node[i].e[6]].connect[0]].loc == 27)
				{
					node[i].loc = 30;
					node[i].ID = 10;
				}
				if(node[CD[node[i].e[1]].connect[2]].loc == 27 && node[CD[node[i].e[1]].connect[0]].loc == 28)
				{
					node[i].loc = 31;
					node[i].ID = 10;
				}
				if(node[CD[node[i].e[2]].connect[1]].loc == 28 && node[CD[node[i].e[2]].connect[3]].loc == 27)
				{
					node[i].loc = 32;
					node[i].ID = 10;
				}	*/
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
				{
					node[i].loc = 33;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 34;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 35;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 36;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 37;
					node[i].corner_ID = 11;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 38;
					node[i].corner_ID = 12;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 39;
					node[i].corner_ID = 13;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 40;
					node[i].corner_ID = 14;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 41;
				}
				if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 42;
				}
				if (node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 43;
				}
				if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 44;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 45;
					node[i].corner_ID = 15;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 46;
					node[i].corner_ID = 16;
				}
				if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 47;
					node[i].corner_ID = 17;
				}
				if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 48;
					node[i].corner_ID = 18;
				}	
				
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
				{
					node[i].loc = 49;
				}
				if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
				{
					node[i].loc = 50;
				}
				if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 51;
				}
				if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
				{
					node[i].loc = 52;
				}				
			}
			/**********************************************************************************************/
			/***********************************node connectivity******************************************/
//			#pragma omp parallel private(j,element, n1) shared(node, CD)
//			{	
//				#pragma omp for
				for(i = 1; i <= NUMNP; i++)
				{
					if (node[i].e[0] != 0)
					{
						node[i].n_n[0] = CD[node[i].e[0]].connect[6];
						node[i].n_n[3] = CD[node[i].e[0]].connect[4];
						node[i].n_n[4] = CD[node[i].e[0]].connect[1];
					}
					if (node[i].e[1] != 0)
					{
						node[i].n_n[0] = CD[node[i].e[1]].connect[2];
						node[i].n_n[3] = CD[node[i].e[1]].connect[0];
						node[i].n_n[5] = CD[node[i].e[1]].connect[5];
					}
					if (node[i].e[2] != 0)
					{
						node[i].n_n[0] = CD[node[i].e[2]].connect[3];
						node[i].n_n[1] = CD[node[i].e[2]].connect[1];
						node[i].n_n[5] = CD[node[i].e[2]].connect[4];
					}
					if (node[i].e[3] != 0)
					{
						node[i].n_n[0] = CD[node[i].e[3]].connect[7];
						node[i].n_n[1] = CD[node[i].e[3]].connect[5];
						node[i].n_n[4] = CD[node[i].e[3]].connect[0];
					}
					if (node[i].e[4] != 0)
					{
						node[i].n_n[2] = CD[node[i].e[4]].connect[5];
						node[i].n_n[3] = CD[node[i].e[4]].connect[7];
						node[i].n_n[4] = CD[node[i].e[4]].connect[2];
					}
					if (node[i].e[5] != 0)
					{
						node[i].n_n[2] = CD[node[i].e[5]].connect[1];
						node[i].n_n[3] = CD[node[i].e[5]].connect[3];
						node[i].n_n[5] = CD[node[i].e[5]].connect[6];
					}
					if (node[i].e[6] != 0)
					{
						node[i].n_n[2] = CD[node[i].e[6]].connect[0];
						node[i].n_n[1] = CD[node[i].e[6]].connect[2];
						node[i].n_n[5] = CD[node[i].e[6]].connect[7];
					}
					if (node[i].e[7] != 0)
					{
						node[i].n_n[2] = CD[node[i].e[7]].connect[4];
						node[i].n_n[1] = CD[node[i].e[7]].connect[6];
						node[i].n_n[4] = CD[node[i].e[7]].connect[3];
					}
				/*	if (node[i].loc >=27 && node[i].loc <=32)
					{
						node[i].n_n[1] = 0;
					} */
				}
				
				/****** special case extended Loc....carefull**********************/
//				#pragma omp for 
				for(i = 1; i <= NUMNP; i++)
				{
					if (node[i].ID == 3 && node[i].loc == 0 && node[i].n_n[3] != 0 && ((node[i].n_n[4] != 0 && node[node[i].n_n[4]].ID == 0) || (node[i].n_n[5] != 0 && node[node[i].n_n[5]].ID == 0)))
					{
						node[i].loc = 27;
						node[i].ID = 10;
						if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
						{
							node[i].n_n[1] = 0;
						}
					}
				/*	if (node[i].ID == 3 && node[i].loc == 0 && node[i].n_n[3] != 0 && ((node[i].n_n[0] != 0 && node[node[i].n_n[0]].ID == 0) || (node[i].n_n[2] != 0 && node[node[i].n_n[2]].ID == 0)))
					{
						node[i].loc = 28;
						node[i].ID = 10;
						if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
						{
							node[i].n_n[1] = 0;
						}
					} */
					if (node[i].ID == 3 && node[i].n_n[3] != 0 && ((node[i].n_n[0] !=0 && node[i].n_n[2] != 0)) && (node[node[i].n_n[4]].ID == 3 || node[node[i].n_n[5]].ID == 3))
					{
						node[i].loc = 28;
						node[i].ID = 10;
						if (node[node[i].n_n[3]].loc == 0)
						{
							node[i].n_n[1] = 0;
						}
					}		
					
					if (node[i].ID == 3 && node[node[i].n_n[5]].loc == 0 && node[i].n_n[5] != 0)  
					{
						node[i].loc = 6;
					}
					
					if (node[i].ID == 3 && node[node[i].n_n[4]].loc == 0 && node[i].n_n[4] != 0)
					{
						node[i].loc = 5;
					}
				}
//			}
			
		/*	for(i=1; i<=NUMNP; i++)
			{
				if (node[i].loc >=27 && node[i].loc <=32)
				{
					node[i].n_n[1] = 0;
				}
			}
			*/
			for(i=1; i<=NUMNP; i++)
			{
			/*	if (node[i].ID == 3 && ((node[node[i].n_n[2]].loc == 27 && (node[i].n_n[5] == 0 || node[i].n_n[5] != 0)) || (node[node[i].n_n[4]].loc == 28 && (node[i].n_n[5] == 0 || node[i].n_n[5] != 0))))
				{
					node[i].loc = 29;
					node[i].ID = 10;
					if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
					{
						node[i].n_n[1] = 0;
					}
				}
				if (node[i].ID == 3 && ((node[node[i].n_n[2]].loc == 27  && (node[i].n_n[4] == 0 || node[i].n_n[4] != 0)) || (node[node[i].n_n[5]].loc == 28 && (node[i].n_n[4] == 0 || node[i].n_n[4] != 0))))
				{
					node[i].loc = 30;
					node[i].ID = 10;
					if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
					{
						node[i].n_n[1] = 0;
					}
				}
				if (node[i].ID == 3 && ((node[node[i].n_n[0]].loc == 27  && (node[i].n_n[4] == 0 || node[i].n_n[4] != 0)) || (node[node[i].n_n[5]].loc == 28 && (node[i].n_n[4] == 0 || node[i].n_n[4] != 0))))
				{
					node[i].loc = 31;
					node[i].ID = 10;
					if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
					{
						node[i].n_n[1] = 0;
					}
				}
				if (node[i].ID == 3 && ((node[node[i].n_n[0]].loc == 27 && (node[i].n_n[5] == 0 || node[i].n_n[5] != 0)) || (node[node[i].n_n[4]].loc == 28 && (node[i].n_n[5] == 0 || node[i].n_n[5] != 0))))
				{
					node[i].loc = 32;
					node[i].ID = 10;
					if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
					{
						node[i].n_n[1] = 0;
					}
				} 
			*/
				if (node[i].loc == 0 && node[node[i].n_n[0]].loc == 1)
				{
					node[i].val = 119;
				}
				if (node[i].loc == 0 && node[node[i].n_n[3]].loc == 4)
				{
					node[i].val = 117;
				}
				if (node[i].loc == 0 && node[node[i].n_n[1]].loc == 2)
				{
					node[i].val = 115;
				}
				if (node[i].loc == 0 && node[node[i].n_n[5]].loc == 5)
				{
					node[i].val = 116;
				}
				if (node[i].loc == 0 && node[node[i].n_n[4]].loc == 6)
				{
					node[i].val = 118;
				}
				if (node[i].loc == 0 && node[node[i].n_n[2]].loc == 3)
				{
					node[i].val = 120;
				}
			}
			
			for(i=1; i<=NUMNP; i++)
			{
				if (node[i].loc == 0)
				{			
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[5]].loc == 5 )
					{ 
						node[i].val =107;
					}
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[1]].loc == 2)
					{
						node[i].val =108;
					}
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[4]].loc == 6)
					{
						node[i].val =109;
					}
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[3]].loc == 4)
					{
						node[i].val =110;
					}
					
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[5]].loc == 5)
					{
						node[i].val =111;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[1]].loc == 2)
					{
						node[i].val =112;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[4]].loc == 6)
					{
						node[i].val =113;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[3]].loc == 4)
					{
						node[i].val =114;
					}
						
					
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[4]].loc == 6 )
					{
						node[i].val=99;
					}
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[5]].loc == 5 )
					{
						node[i].val=100;
					}
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[5]].loc == 5 )
					{
						node[i].val=101;
					}
					if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[4]].loc == 6 )
					{
						node[i].val=102;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[4]].loc == 6 )
					{
						node[i].val=103;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[5]].loc == 5 )
					{
						node[i].val=104;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[5]].loc == 5 )
					{
						node[i].val=105;
					}
					if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[4]].loc == 6 )
					{
						node[i].val=106;
					}
				}
				
				if (node[i].loc != 0)
				{	
					if (node[node[i].n_n[2]].loc == 15 && node[node[i].n_n[3]].loc == 23 )
					{ 
						node[i].val =43;
					}
					if (node[node[i].n_n[2]].loc == 15 && node[node[i].n_n[1]].loc == 26)
					{
						node[i].val =50;
					}
					if (node[node[i].n_n[1]].loc == 26 && node[node[i].n_n[0]].loc == 19)
					{
						node[i].val =66;
					}
					if (node[node[i].n_n[0]].loc == 19 && node[node[i].n_n[3]].loc == 23)
					{
						node[i].val =59;
					}
						
					if (node[node[i].n_n[2]].loc == 16 && node[node[i].n_n[4]].loc == 23)
					{
						node[i].val =44;
					}
					if (node[node[i].n_n[4]].loc == 23 && node[node[i].n_n[0]].loc == 20)
					{
						node[i].val =60;
					}
					if (node[node[i].n_n[5]].loc == 24 && node[node[i].n_n[0]].loc == 20)
					{
						node[i].val =61;
					}
					if (node[node[i].n_n[5]].loc == 24 && node[node[i].n_n[2]].loc == 16)
					{
						node[i].val =45;
					}
					
					if (node[node[i].n_n[2]].loc == 17 && node[node[i].n_n[3]].loc == 24)
					{
						node[i].val =46;
					}
					if (node[node[i].n_n[3]].loc == 24 && node[node[i].n_n[0]].loc == 21)
					{
						node[i].val =62;
					}
					if (node[node[i].n_n[0]].loc == 21 && node[node[i].n_n[1]].loc == 25)
					{
						node[i].val =63;
					}
					if (node[node[i].n_n[2]].loc == 17 && node[node[i].n_n[1]].loc == 25)
					{
						node[i].val =47;
					}
						
					if (node[node[i].n_n[5]].loc == 25 && node[node[i].n_n[0]].loc == 22)
					{
						node[i].val =64;
					}
					if (node[node[i].n_n[4]].loc == 26 && node[node[i].n_n[0]].loc == 22)
					{
						node[i].val =65;
					}
					if (node[node[i].n_n[4]].loc == 26 && node[node[i].n_n[2]].loc == 18)
					{
						node[i].val =49;
					}
					if (node[node[i].n_n[5]].loc == 25 && node[node[i].n_n[2]].loc == 18)
					{
						node[i].val =48;
					}
						
					if (node[node[i].n_n[4]].loc == 19 && node[node[i].n_n[3]].loc == 20)
					{
						node[i].val =67;
					}
					if (node[node[i].n_n[3]].loc == 20 && node[node[i].n_n[5]].loc == 21)
					{
						node[i].val =68;
					}
					if (node[node[i].n_n[5]].loc == 21 && node[node[i].n_n[1]].loc == 22)
					{
						node[i].val =69;
					}
					if (node[node[i].n_n[1]].loc == 22 && node[node[i].n_n[4]].loc == 19)
					{
						node[i].val =70;
					}
					
					if (node[node[i].n_n[4]].loc == 15 && node[node[i].n_n[3]].loc == 16)
					{
						node[i].val =71;
					}
					if (node[node[i].n_n[3]].loc == 16 && node[node[i].n_n[5]].loc == 17)
					{
						node[i].val =72;
					}
					if (node[node[i].n_n[5]].loc == 17 && node[node[i].n_n[1]].loc == 18)
					{
						node[i].val =73;
					}
					if (node[node[i].n_n[4]].loc == 15 && node[node[i].n_n[1]].loc == 18)
					{
						node[i].val =74;
					}
					
					
					if (node[node[i].n_n[2]].loc == 16 && node[i].loc == 4 && node[node[i].n_n[4]].loc == 4 && node[node[i].n_n[5]].loc == 4)
					{
						node[i].val = 86;
					}
					if (node[node[i].n_n[2]].loc == 16 && node[i].loc == 3 && node[node[i].n_n[4]].loc == 3 && node[node[i].n_n[5]].loc == 3)
					{
						node[i].val = 77;
					}
		
					if (node[node[i].n_n[2]].loc == 17 && node[i].loc == 26 && node[node[i].n_n[1]].loc == 5 && node[node[i].n_n[3]].loc == 5)
					{
						node[i].val = 85;
					}
					if (node[node[i].n_n[5]].loc == 17 && node[i].loc == 3 && node[node[i].n_n[1]].loc == 3 && node[node[i].n_n[3]].loc == 3)
					{
						node[i].val = 93;
					}

					if (node[node[i].n_n[2]].loc == 18 && node[i].loc == 2 && node[node[i].n_n[4]].loc == 2 && node[node[i].n_n[5]].loc == 2)
					{
						node[i].val = 84;
					}					
					if (node[node[i].n_n[1]].loc == 18 && node[i].loc == 3 && node[node[i].n_n[4]].loc == 3 && node[node[i].n_n[5]].loc == 3)
					{
						node[i].val = 81;
					}					

					if (node[node[i].n_n[4]].loc == 19 && node[i].loc == 1 && node[node[i].n_n[1]].loc == 1 && node[node[i].n_n[3]].loc == 1)
					{
						node[i].val = 95;
					}
					if (node[node[i].n_n[0]].loc == 19 && node[i].loc == 6 && node[node[i].n_n[1]].loc == 6 && node[node[i].n_n[3]].loc == 6)
					{
						node[i].val = 87;
					}

					if (node[node[i].n_n[3]].loc == 20 && node[i].loc == 1 && node[node[i].n_n[4]].loc == 1 && node[node[i].n_n[5]].loc == 1)
					{
						node[i].val = 75;
					}
					if (node[node[i].n_n[0]].loc == 20 && node[i].loc == 4 && node[node[i].n_n[4]].loc == 4 && node[node[i].n_n[5]].loc == 4)
					{
						node[i].val = 90;
					}

					if (node[node[i].n_n[5]].loc == 21 && node[i].loc == 1 && node[node[i].n_n[1]].loc == 1 && node[node[i].n_n[3]].loc == 1)
					{
						node[i].val = 91;
					}
					if (node[node[i].n_n[0]].loc == 21 && node[i].loc == 5 && node[node[i].n_n[1]].loc == 5 && node[node[i].n_n[3]].loc == 5)
					{
						node[i].val = 89;
					}

					if (node[node[i].n_n[1]].loc == 22 && node[i].loc == 1 && node[node[i].n_n[4]].loc == 1 && node[node[i].n_n[5]].loc == 1)
					{
						node[i].val = 79;
					}
					if (node[node[i].n_n[0]].loc == 22 && node[i].loc == 3 && node[node[i].n_n[4]].loc == 3 && node[node[i].n_n[5]].loc == 3)
					{
						node[i].val = 88;
					}

					if (node[node[i].n_n[4]].loc == 23 && node[i].loc == 4 && node[node[i].n_n[0]].loc == 4 && node[node[i].n_n[2]].loc == 4)
					{
						node[i].val = 96;
					}
					if (node[node[i].n_n[3]].loc == 23 && node[i].loc == 6 && node[node[i].n_n[0]].loc == 6 && node[node[i].n_n[2]].loc == 6)
					{
						node[i].val = 78;
					}

					if (node[node[i].n_n[5]].loc == 24 && node[i].loc == 4 && node[node[i].n_n[0]].loc == 4 && node[node[i].n_n[2]].loc == 4)
					{
						node[i].val = 94;
					}
					if (node[node[i].n_n[3]].loc == 24 && node[i].loc == 5 && node[node[i].n_n[0]].loc == 5 && node[node[i].n_n[2]].loc == 5)
					{
						node[i].val = 76;
					}

					if (node[node[i].n_n[1]].loc == 25 && node[i].loc == 5 && node[node[i].n_n[0]].loc == 5 && node[node[i].n_n[2]].loc == 5)
					{
						node[i].val = 80;
					}
					if (node[node[i].n_n[5]].loc == 25 && node[i].loc == 2 && node[node[i].n_n[0]].loc == 2 && node[node[i].n_n[2]].loc == 2)
					{
						node[i].val = 92;
					}

					if (node[node[i].n_n[4]].loc == 26 && node[i].loc == 2 && node[node[i].n_n[0]].loc == 2 && node[node[i].n_n[2]].loc == 2)
					{
						node[i].val = 98;
					}
					if (node[node[i].n_n[1]].loc == 26 && node[i].loc == 6 && node[node[i].n_n[0]].loc == 6 && node[node[i].n_n[2]].loc == 6)
					{
						node[i].val = 82;
					}	
						
					if (node[node[i].n_n[3]].loc == 12 && node[i].loc == 17)
					{
						node[i].val = 38;
					}
					if (node[node[i].n_n[5]].loc == 12 && node[i].loc == 16)
					{
						node[i].val = 37;
					}
					
					if (node[node[i].n_n[4]].loc == 11 && node[i].loc == 16)
					{
						node[i].val = 36;
					}
					if (node[node[i].n_n[3]].loc == 11 && node[i].loc == 15)
					{
						node[i].val = 35;
					}				
						
					if (node[node[i].n_n[1]].loc == 14 && node[i].loc == 15)
					{
						node[i].val = 42;
					}
					if (node[node[i].n_n[4]].loc == 14 && node[i].loc == 18)
					{
						node[i].val = 41;
					}
						
					if (node[node[i].n_n[5]].loc == 13 && node[i].loc == 18)
					{
						node[i].val = 40;
					}
					if (node[node[i].n_n[3]].loc == 13 && node[i].loc == 17)
					{
						node[i].val = 39;
					}
						
					if (node[node[i].n_n[3]].loc == 8 && node[i].loc == 21)
					{
						node[i].val = 54;
					}
					if (node[node[i].n_n[5]].loc == 8 && node[i].loc == 20)
					{
						node[i].val = 53;
					}
						
					if (node[node[i].n_n[4]].loc == 7 && node[i].loc == 20)
					{
						node[i].val = 52;
					}
					if (node[node[i].n_n[3]].loc == 7 && node[i].loc == 19)
					{
						node[i].val = 51;
					}
						
					if (node[node[i].n_n[1]].loc == 10 && node[i].loc == 19)
					{
						node[i].val = 58;
					}
					if (node[node[i].n_n[4]].loc == 10 && node[i].loc == 22)
					{
						node[i].val = 57;
					}
					
					if (node[node[i].n_n[1]].loc == 9 && node[i].loc == 21)
					{
						node[i].val = 55;
					}
					if (node[node[i].n_n[5]].loc == 9 && node[i].loc == 22)
					{
						node[i].val = 56;
					}
						
					if (node[i].n_n[3] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[0]].loc == 7)
					{
						node[i].val = 27;
					}
					if (node[i].n_n[3] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[0]].loc == 8)
					{
						node[i].val = 28;
					}
					if (node[i].n_n[1] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[0]].loc == 9)
					{
						node[i].val = 29;
					}
					if (node[i].n_n[1] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[0]].loc == 10)
					{
						node[i].val = 30;
					}
						
					if (node[i].n_n[3] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[2]].loc == 11)
					{
						node[i].val = 31;
					}
					if (node[i].n_n[3] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[2]].loc == 12)
					{
						node[i].val = 32;
					}
					if (node[i].n_n[1] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[2]].loc == 13)
					{
						node[i].val = 33;
					}
					if (node[i].n_n[1] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[2]].loc == 14)
					{
						node[i].val = 34;
					}
				}												
			}
			/***** val near singular points **************/
			for (i=1; i<= NUMNP; i++)
			{			
				if (node[i].loc == 27 || node[i].loc == 28 || node[i].loc == 29 || node[i].loc == 30 || node[i].loc == 31 || node[i].loc == 32) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[4]].val = 116;
				}
				if (node[i].loc == 27 ) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[1]].val = 117;
					node[node[i].n_n[3]].val = 115;
				}
				if (node[i].loc == 28 ) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[0]].val = 120;
					node[node[i].n_n[2]].val = 119;
				}
				if (node[i].loc == 29 ) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[0]].val = 120;
					node[node[i].n_n[4]].val = 116;
					node[node[i].n_n[1]].val = 115;
				}
				if (node[i].loc == 30 ) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[0]].val = 120;
					node[node[i].n_n[1]].val = 117;
					node[node[i].n_n[4]].val = 116;
				}
				if (node[i].loc == 31 ) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[2]].val = 119;
					node[node[i].n_n[1]].val = 117;
					node[node[i].n_n[4]].val = 116;
				}
				if (node[i].loc == 32 ) 
				{
					node[i].corner_ID = 0;
					node[node[i].n_n[2]].val = 119;
					node[node[i].n_n[3]].val = 115;
					node[node[i].n_n[4]].val = 116;
				}
				if (node[node[i].n_n[4]].loc == 27 && node[i].n_n[1] == 0)
				{
					node[i].val = 96;
				}
				if (node[node[i].n_n[4]].loc == 27 && node[i].n_n[3] == 0)
				{
					node[i].val = 98;
				}
				if (node[node[i].n_n[4]].loc == 28 && node[i].n_n[0] == 0)
				{
					node[i].val = 95;
				}
				if (node[node[i].n_n[4]].loc == 28 && node[i].n_n[2] == 0)
				{
					node[i].val = 97;
				}
				
				if (node[node[i].n_n[4]].loc == 29 && node[i].loc == 20)
				{
					node[i].val = 52;
				}
				if (node[node[i].n_n[4]].loc == 30 && node[i].loc == 22)
				{
					node[i].val = 57;
				}
				if (node[node[i].n_n[4]].loc == 31 && node[i].loc == 18)
				{
					node[i].val = 41;
				}
				if (node[node[i].n_n[4]].loc == 32 && node[i].loc == 16)
				{
					node[i].val = 36;
				}
				
				if (node[node[i].n_n[4]].loc == 29 && node[i].loc == 18)
				{
					node[i].val = 41;
				}
				if (node[node[i].n_n[4]].loc == 30 && node[i].loc == 16)
				{
					node[i].val = 36;
				}
				if (node[node[i].n_n[4]].loc == 31 && node[i].loc == 20)
				{
					node[i].val = 52;
				}
				if (node[node[i].n_n[4]].loc == 32 && node[i].loc == 22)
				{
					node[i].val = 57;
				}
				
				if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[0]].loc == 20)
				{
					node[i].val = 60;
				}
				if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[0]].loc == 22)
				{
					node[i].val = 65;
				}
				if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[2]].loc == 16)
				{
					node[i].val = 44;
				}
				if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[2]].loc == 18)
				{
					node[i].val = 49;
				}
				
				if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[3]].loc == 20)
				{
					node[i].val = 67;
				}
				if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[1]].loc == 22)
				{
					node[i].val = 70;
				}
				if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[3]].loc == 16)
				{
					node[i].val = 71;
				}
				if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[1]].loc == 18)
				{
					node[i].val = 74;
				}
				if (node[node[i].n_n[3]].loc == 27 && node[node[i].n_n[0]].loc == 28)
				{
					node[i].val = 59;
				}
				if (node[node[i].n_n[1]].loc == 27 && node[node[i].n_n[0]].loc == 28)
				{
					node[i].val = 66;
				}
				if (node[node[i].n_n[3]].loc == 27 && node[node[i].n_n[2]].loc == 28)
				{
					node[i].val = 43;
				}
				if (node[node[i].n_n[1]].loc == 27 && node[node[i].n_n[2]].loc == 28)
				{
					node[i].val = 50;
				}
				
				if (node[node[i].n_n[3]].loc == 29 && node[i].loc == 19)
				{
					node[i].val = 51;
				}
				if (node[node[i].n_n[1]].loc == 30 && node[i].loc == 19)
				{
					node[i].val = 58;
				}
				if (node[node[i].n_n[3]].loc == 32 && node[i].loc == 15)
				{
					node[i].val = 35;
				}
				if (node[node[i].n_n[1]].loc == 31 && node[i].loc == 15)
				{
					node[i].val = 42;
				}
				if (node[node[i].n_n[0]].loc == 29)
				{
					node[i].val = 27;
				}
				if (node[node[i].n_n[0]].loc == 30)
				{
					node[i].val = 30;
				}
				if (node[node[i].n_n[2]].loc == 31)
				{
					node[i].val = 34;
				}
				if (node[node[i].n_n[2]].loc == 32)
				{
					node[i].val = 31;
				}
			}
			
			all_bou_node = inl+wal+bou+out+corn_node;
			
			for (i=1; i<= NUMNP; i++)
			{
				if (node[i].ID == 10)
				{
					node[i].n_n[1] = 0;
				}
				
			}
			
			/********************************writing node neighbour file*****************************************/
			printf("creating node_neighbour.neu\n");
	/*		fnode = fopen("node_neighbour.neu","w");
			for(i=1; i<=NUMNP; i++)
			{
				fprintf(fnode,"%e %e %e %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", node[i].x, node[i].y, node[i].z, node[i].n_n[0], node[i].n_n[1], node[i].n_n[2], node[i].n_n[3], node[i].n_n[4], node[i].n_n[5], node[i].e[0], node[i].e[1], node[i].e[2], node[i].e[3], node[i].e[4], node[i].e[5], node[i].e[6], node[i].e[7], node[i].loc, node[i].val, node[i].ID, node[i].corner_ID);
			}		
			fclose(fnode);	
	*/

			position=0;
			MPI_Pack_size((NUMNP+10)*3,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
			Bcast_buffer1 = malloc(memsize);
			for(i=1; i<=NUMNP; i++)
			{
				MPI_Pack(&node[i].x, 1, MPI_DOUBLE, Bcast_buffer1, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].y, 1, MPI_DOUBLE, Bcast_buffer1, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].z, 1, MPI_DOUBLE, Bcast_buffer1, memsize, &position, MPI_COMM_WORLD);
				
			}	
			
			position=0;
			MPI_Pack_size((NUMNP+10)*10,MPI_INT,MPI_COMM_WORLD,&memsize);
			Bcast_buffer2 = malloc(memsize);
			for(i=1; i<=NUMNP; i++)
			{
				MPI_Pack(&node[i].n_n[0], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].n_n[1], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].n_n[2], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].n_n[3], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].n_n[4], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].n_n[5], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].loc, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].ID, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].proc, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
				MPI_Pack(&node[i].corner_ID, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			}	
			printf("done broadcasting\n");
			
			temp = 0;
			fnode = fopen("boundary_count.neu","w");
			for(i=0; i<=all_bou_node; i++)
			{
				if (i == 1 || i == inl || i == inl+out || i == inl+out+wal || i == inl+out+wal+bou)
				{
					temp = 0;
				}
				if (i == 0)
				{
					fprintf(fnode,"%d %d %d %d %d\n", inl, out, bou, wal, corn_node);
				}
				if (i <= inl && i != 0)
				{
					fprintf(fnode,"%d\n",inlet_node[temp]);
					temp++;
				}
				if (i <= inl+out && i > inl)
				{
					fprintf(fnode,"%d\n",outlet_node[temp]);
					temp++;
				}
				if (i <= inl+out+wal && i > inl+out)
				{
					fprintf(fnode,"%d\n",wall_node[temp]);
					temp++;
				}
				if (i <= inl+out+wal+bou && i > inl+out+wal)
				{
					fprintf(fnode,"%d\n",boundary_node[temp]);
					temp++;
				}
				if (i <= inl+out+wal+bou+corn_node && i > inl+out+wal+bou)
				{
					fprintf(fnode,"%d\n",corner_node[temp]);
					temp++;
				}
			}		
			fclose(fnode);
			
			fnode = fopen("metis_input.txt","w");
			fprintf(fnode,"%d 3\n",NELEM);
			for(i=1; i<=NELEM; i++)
			{
				fprintf(fnode,"%d %d %d %d %d %d %d %d\n", CD[i].connect[4], CD[i].connect[0], CD[i].connect[1], CD[i].connect[5] , CD[i].connect[7], CD[i].connect[3], CD[i].connect[2], CD[i].connect[6]);
			}		
			fclose(fnode);
		
			printf("started domain decomposition using metis\n");
		//	system("cp metis_input.txt $HOME/metis-4.0/.");
			sprintf(filename,"./partdmesh metis_input.txt %d",size-1);
			system(filename);
			sprintf(filename,"mv metis_input.txt.npart.%d nodes_proc.txt",size-1);
			system(filename);
			printf("completed decomposition using metis\n");
		
			felem = fopen("elem_neighbour.neu","w");
			for(i=1; i<=NELEM; i++)
			{
				fprintf(felem,"%d %d %d %d %d %d %d %d\n",CD[i].connect[0], CD[i].connect[1], CD[i].connect[2], CD[i].connect[3], CD[i].connect[4], CD[i].connect[5], CD[i].connect[6], CD[i].connect[7]);
			}		
			fclose(felem);
		
			fnode = fopen("nodes_proc.txt","rt");
			for(i=1; i<= NUMNP; i++)
			{
				fscanf(fnode,"%d\n",&node[i].proc);
				node[i].proc = node[i].proc+1;
			}
			fclose(fnode);
		
			for (i=1;i<size;i++)
			{
				sprintf(filename,"nodes_proc_%d.txt",i);
				fnode=fopen(filename,"w");
				for (j=1;j<=NUMNP;j++)
				{
					fprintf(fnode,"%d\n",node[j].proc);
				}
				fclose(fnode);
			}
		
			/*************************/
			sprintf(filename,"nodefile.plt");
			fp= fopen(filename,"w");
			fprintf(fp,"TITLE = \"Node file\"\n");
			fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"Proc\", \"Loc\", \n");
			fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",NUMNP,NELEM);
			for(i=0;i<NUMNP;i++)
			{
				fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\n", node[i+1].x, node[i+1].y, node[i+1].z, node[i+1].proc, node[i+1].loc);
			}
			fprintf(fp,"\n\n\n");
	
			for(i=1; i<=NELEM; i++)
			{
				fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CD[i].connect[0], CD[i].connect[1], CD[i].connect[5], CD[i].connect[4], CD[i].connect[3], CD[i].connect[2], CD[i].connect[6], CD[i].connect[7] );
			}
			fclose(fp);
			/********************/
		}
		/****************************************Reading neighbour file***********************************************/
		printf("\n\n\n");	

		/**********************************************************************************************
				  storing all boundary nodes in a single array.... 
	              this array will be freed before solver is executed
		***********************************************************************************************/			
		BOUND *temp_check;
		bound=0;
		all_bou_node = inl+out+wal+bou;
		all_boundary_nodes=(int*)malloc((all_bou_node+10)*sizeof(int));
		temp_check = (BOUND *)malloc((NUMNP+10)*sizeof(BOUND));
		memset(all_boundary_nodes, 0, (all_bou_node+10)*sizeof(int));
		memset(temp_check, 0, (NUMNP+10)*sizeof(int));
		
		for (i=0;i<inl;i++)
		{
			if (node[inlet_node[i]].corner_ID ==0 && temp_check[inlet_node[i]].check == 0)
			{ 
				all_boundary_nodes[bound]=inlet_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i=0;i<out;i++)
		{
			if (node[outlet_node[i]].corner_ID ==0 && temp_check[outlet_node[i]].check == 0)
			{ 
				all_boundary_nodes[bound]=outlet_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i=0;i<wal;i++)
		{
			if (node[wall_node[i]].corner_ID ==0 && temp_check[wall_node[i]].check == 0)
			{ 
				all_boundary_nodes[bound]=wall_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i=0;i<bou;i++)
		{
			if (node[boundary_node[i]].corner_ID ==0 && temp_check[boundary_node[i]].check == 0)
			{ 
				all_boundary_nodes[bound] = boundary_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i = 0; i < corn_node; i++)
		{
			if (temp_check[corner_node[i]].check == 0)
			{ 
				all_boundary_nodes[bound] = corner_node[i];
				bound++;
			}
		}
		all_bou_node = bound;
		/***********************jacobian and Metrics calculation**************************************/
		//all_bou=bound;
		memory=(NUMNP+1)+(3*all_bou_node);
		memory_elem=((NELEM+1)+(3*all_bou_node));
		metric = (TRANSFORMED *)malloc(memory*sizeof(TRANSFORMED));
		jacobian = (JACOB *)malloc(memory*sizeof(JACOB));
		deter = (DETERM *)malloc(memory*sizeof(DETERM));
		det = (double *)malloc(memory*sizeof(double));
		
		jac = (JAC *)malloc(memory*sizeof(JAC));
		jaco = (TEMP_JAC *)malloc(memory*sizeof(TEMP_JAC));
		
	/*	memset(jacobian, 0, sizeof(jacobian));
		memset(metric, 0, sizeof(metric));
		memset(det,0, sizeof(det));
		*/
		if (jacobian == NULL)
		{
			printf("memory not allocatted\n");
		}
	
		/***************************************************************************************************/		
		metric_term(CD, node, NELEM, jacobian, metric, cor, NUMNP, det, outlet_node, out_node, boundary_node, bou_node, inlet_node, inl_node, wall_node, wal_node, all_boundary_nodes, all_bou_node, jac, jaco);
		
		free(jac);
		jac = NULL;
		free(jaco);
		jaco = NULL;
		free(temp_check);
		temp_check = NULL;
		/*************************/
		sprintf(filename,"metric_loc.plt");
		fp= fopen(filename,"w");
		fprintf(fp,"TITLE = \"Node file\"\n");
		fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"Proc\", \"val\", \"loc\", \"det\",\n");
		fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",NUMNP,NELEM);
		for(i=0;i<NUMNP;i++)
		{
			fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf\n", node[i+1].x, node[i+1].y, node[i+1].z, node[i+1].proc, node[i+1].ID, node[i+1].loc, 1.0/det[i+1]);
		}
		fprintf(fp,"\n\n\n");
	
		for(i=1; i<=NELEM; i++)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CD[i].connect[0], CD[i].connect[1], CD[i].connect[5], CD[i].connect[4], CD[i].connect[3], CD[i].connect[2], CD[i].connect[6], CD[i].connect[7] );
		}
		fclose(fp);
		/********************/
		
		
		/***************SEGREGATING AND WRITING BOUNDARY NODES FOR CORRESPONDING PROCESSOR******************/
		
		for (i=1; i<size; i++)
		{
			inl =0;
			out =0;
			bou =0;
			wal =0;
			for(j=0;j<inl_node;j++)
			{
				if (node[inlet_node[j]].proc == i)
				{
					inl++;
				}
			}
			sprintf(filename,"inlet_%d.txt",i);
			fnode = fopen(filename,"w");
			for(j=0;j<inl_node;j++)
			{
				if (node[inlet_node[j]].proc == i)
				{
					fprintf(fnode,"%d\n",inlet_node[j]);
					
				}
			}
			fclose(fnode);

			for(j=0;j<out_node;j++)
			{
				if (node[outlet_node[j]].proc == i)
				{
					out++;
				}
			}
			sprintf(filename,"outlet_%d.txt",i);
			fnode = fopen(filename,"w");
			for(j=0;j<out_node;j++)
			{
				if (node[outlet_node[j]].proc == i)
				{
					fprintf(fnode,"%d\n",outlet_node[j]);
				}
			}
			fclose(fnode);

			for(j=0;j<wal_node;j++)
			{
				if (node[wall_node[j]].proc == i)
				{
					wal++;
				}
			}
			sprintf(filename,"wall_%d.txt",i);
			fnode = fopen(filename,"w");
			for(j=0;j<wal_node;j++)
			{
				if (node[wall_node[j]].proc == i)
				{
					fprintf(fnode,"%d\n",wall_node[j]);
				}
			}
			fclose(fnode);

			for(j=0;j<bou_node;j++)
			{
				if (node[boundary_node[j]].proc == i)
				{
					bou++;
				}
			}
			sprintf(filename,"boundary_%d.txt",i);
			fnode = fopen(filename,"w");
			for(j=0;j<bou_node;j++)
			{
				if (node[boundary_node[j]].proc == i)
				{
					fprintf(fnode,"%d\n",boundary_node[j]);
				}
			}
			fclose(fnode);
			
			MPI_Send(&inl, 1, MPI_INT, i,i, MPI_COMM_WORLD);
			MPI_Send(&out, 1, MPI_INT, i,i, MPI_COMM_WORLD);
			MPI_Send(&wal, 1, MPI_INT, i,i, MPI_COMM_WORLD);
			MPI_Send(&bou, 1, MPI_INT, i,i, MPI_COMM_WORLD);
		}
		/*******************************give boundary ID to nodes******************************************/
		for (i=0; i<out_node; i++)
		{
			node[outlet_node[i]].ID = 2;
		}
		for (i=0; i<inl_node; i++)
		{
			node[inlet_node[i]].ID = 1;
		}
		for (i=0; i<wal_node; i++)
		{
			if (node[i].loc != 28)
			{
				node[wall_node[i]].ID = 3;
			}
		}
		for (i=0; i<bou_node; i++)
		{
			node[boundary_node[i]].ID = 4;
		}
		/********************************Write node_location.txt file**************************************/
	/*	fnode = fopen("node_location.txt","w");
		for(i=1; i<=NUMNP; i++)
		{
			fprintf(fnode,"%d %d %d %d\n", node[i].val, node[i].loc, node[i].corner_ID, node[i].ID);
		}		
		fclose(fnode);
		*/
		/*********************************WRITE JACOBIAN DATA TO FILE***************************************/
/*		fnode = fopen("transformation_data.txt","w");
		for(i=1; i<=NUMNP; i++)
		{
			fprintf(fnode,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", jacobian[i].x_zeta, jacobian[i].x_eta, jacobian[i].x_xi, jacobian[i].y_zeta, jacobian[i].y_eta, jacobian[i].y_xi, jacobian[i].z_zeta, jacobian[i].z_eta, jacobian[i].z_xi, metric[i].zeta_x, metric[i].eta_x, metric[i].xi_x, metric[i].zeta_y, metric[i].eta_y, metric[i].xi_y, metric[i].zeta_z, metric[i].eta_z, metric[i].xi_z, det[i]);
		}		
		fclose(fnode);
*/		
	
		MPI_Pack_size((NUMNP+10)*5, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		Bcast_buffer3 = malloc(memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Pack(&jacobian[i].x_zeta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].x_eta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].x_xi, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_zeta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_eta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);			
		}
		
		MPI_Pack_size((NUMNP+10)*4, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		Bcast_buffer4 = malloc(memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Pack(&jacobian[i].y_xi, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_zeta, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_eta, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_xi, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			
		}
		/****************************************************************************************************/
	//	iterations = 100000;
		restart = 2;	
		
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	/***********************************Broadcast NUMNP and NELEM**********************************/
	MPI_Bcast(&NUMNP,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&NELEM,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&iterations, 1, MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&restart, 1, MPI_INT,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Pack_size((NUMNP+10)*3,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
	if (myrank >0)
	{
		Bcast_buffer1 = malloc(memsize);
	}	
	MPI_Bcast(Bcast_buffer1,(NUMNP+10)*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if (myrank == 0)
	{
		free(Bcast_buffer1);
		Bcast_buffer1 = NULL;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Pack_size((NUMNP+10)*10,MPI_INT,MPI_COMM_WORLD,&memsize);
	if (myrank >0)
	{
		Bcast_buffer2 = malloc(memsize);
	}
	MPI_Bcast(Bcast_buffer2,(NUMNP+10)*10,MPI_INT,0,MPI_COMM_WORLD);
	if (myrank == 0)
	{
		free(Bcast_buffer2);
		Bcast_buffer2 = NULL;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Pack_size((NUMNP+10)*5,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
	if (myrank >0)
	{
		Bcast_buffer3 = malloc(memsize);
	}
	MPI_Bcast(Bcast_buffer3,(NUMNP+10)*5,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if (myrank == 0)
	{
		free(Bcast_buffer3);
		Bcast_buffer3 = NULL;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Pack_size((NUMNP+10)*4,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
	if (myrank >0)
	{
		Bcast_buffer4 = malloc(memsize);
	}
	MPI_Bcast(Bcast_buffer4,(NUMNP+10)*4,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if (myrank == 0)
	{
		free(Bcast_buffer4);
		Bcast_buffer4 = NULL;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (myrank == 0)
	{
		printf("done proper broadcasting\n");
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/********************************************************************************************/
	/********************************************************************************************/
	if (myrank>0)
	{
		int sd_gh, sd_node, j, h, gh, d[size+1], temp, *proc_node, temp1, no_of_tip_send, no_of_tip_recv;
		MNODE *d_node, *tmp_node;
		SUB_DOM *sd_gh_node;
		int *c, *tmp_c, **b, **gb, **loc_dat, **recv_b, neigh_pro, glob, loca, *recv_c;
		int sd_inl_node, sd_wal_node, sd_out_node, sd_bou_node, *sd_inlet_node, *sd_outlet_node, *sd_wall_node, *sd_boundary_node;
		char line[1100];
		JACOB *t_jacobian;
		TRANSFORMED *t_metric;
		DOM_TIP *tip;
		DOM_TIP *tip_recv;
		double *t_det;
		
		tip = (DOM_TIP*)malloc(((5000))*sizeof(DOM_TIP));
		tip_recv = (DOM_TIP*)malloc(((5000))*sizeof(DOM_TIP));
		
		FILE *fnode;
		
		MPI_Recv(&inl_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		MPI_Recv(&out_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		MPI_Recv(&wal_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		MPI_Recv(&bou_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		
		d_node = (MNODE*)malloc(((NUMNP+10))*sizeof(MNODE));
		memset(d_node, 0, sizeof(d_node));
		for(i=0; i<=NUMNP; i++)
		{
			d_node[i].x=0;
			d_node[i].y=0;
			d_node[i].z=0;
			d_node[i].e[0]=0;
			d_node[i].e[1]=0;
			d_node[i].e[2]=0;
			d_node[i].e[3]=0;
			d_node[i].e[4]=0;
			d_node[i].e[5]=0;
			d_node[i].e[6]=0;
			d_node[i].e[7]=0;
			d_node[i].ID=0;
			d_node[i].n_n[0]=0;
			d_node[i].n_n[1]=0;
			d_node[i].n_n[2]=0;
			d_node[i].n_n[3]=0;
			d_node[i].n_n[4]=0;
			d_node[i].n_n[5]=0;
		//	d_node[i].location=0;
		//	d_node[i].corner=0;
			d_node[i].corner_ID=0;
			d_node[i].val=0;
			d_node[i].loc=0;
		//	d_node[i].weno=0;
			d_node[i].proc=0;
			d_node[i].local=0;
			d_node[i].global=0;
			d_node[i].req=0;
		}
		MPI_Pack_size((NUMNP+10)*3,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer1, memsize, &position, &d_node[i].x, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer1, memsize, &position, &d_node[i].y, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer1, memsize, &position, &d_node[i].z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		}
		free(Bcast_buffer1);
		Bcast_buffer1 = NULL;
		
		MPI_Pack_size((NUMNP+10)*10,MPI_INT,MPI_COMM_WORLD,&memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[0], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[1], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[2], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[3], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[4], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[5], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].loc, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].ID, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].proc, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].corner_ID, 1, MPI_INT, MPI_COMM_WORLD);
		}
		free(Bcast_buffer2);
		Bcast_buffer2 = NULL;
		
		/*
		fnode = fopen("node_neighbour.neu","rt");		
		i=1;	
		while(fgets(line, 200, fnode) != NULL)
		{
			if(i > NUMNP)
			{
				break;
			}
			sscanf(line,"%lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &d_node[i].x, &d_node[i].y, &d_node[i].z, &d_node[i].n_n[0], &d_node[i].n_n[1], &d_node[i].n_n[2], &d_node[i].n_n[3], &d_node[i].n_n[4], &d_node[i].n_n[5], &d_node[i].e[0], &d_node[i].e[1], &d_node[i].e[2], &d_node[i].e[3], &d_node[i].e[4], &d_node[i].e[5], &d_node[i].e[6], &d_node[i].e[7], &d_node[i].loc, &d_node[i].val, &d_node[i].ID, &d_node[i].corner_ID);
			i++;
		}	
		fclose(fnode);
		*/
		
		/******************************************************************************************************/
		/*******************READING PROC LIST AND COUNTING NODES IN SUBDOMAIN**********************************/
		sd_node=0;
		sprintf(filename,"nodes_proc_%d.txt",myrank);
		fnode=fopen(filename,"rt");
		for (i=1;i<=NUMNP;i++)
		{
			fscanf(fnode,"%d\n",&d_node[i].proc);
			if(d_node[i].proc == myrank)
			{
				sd_node++;				
			}
		}
		fclose(fnode);
		sd_node =sd_node;
		tmp_node = (MNODE*)malloc((sd_node+10)*sizeof(MNODE));
		memset(tmp_node, 0, sizeof(tmp_node));
		
		for(i=0; i<=sd_node+8; i++)
		{
			tmp_node[i].x=0;
			tmp_node[i].y=0;
			tmp_node[i].z=0;
			tmp_node[i].e[0]=0;
			tmp_node[i].e[1]=0;
			tmp_node[i].e[2]=0;
			tmp_node[i].e[3]=0;	
			tmp_node[i].e[4]=0;
			tmp_node[i].e[5]=0;
			tmp_node[i].e[6]=0;
			tmp_node[i].e[7]=0;	
			tmp_node[i].ID=0;
			tmp_node[i].n_n[0]=0;
			tmp_node[i].n_n[1]=0;
			tmp_node[i].n_n[2]=0;
			tmp_node[i].n_n[3]=0;
			tmp_node[i].n_n[4]=0;
			tmp_node[i].n_n[5]=0;
		//	tmp_node[i].location=0;
		//	tmp_node[i].corner=0;
			tmp_node[i].corner_ID=0;
			tmp_node[i].val=0;
			tmp_node[i].loc=0;
		//	tmp_node[i].weno=0;
			tmp_node[i].proc=0;
			tmp_node[i].local=0;
			tmp_node[i].global=0;
			tmp_node[i].req=0;
		}
		
		j = 1;
		for (i=1; i<=NUMNP; i++)
		{
			if (d_node[i].proc == myrank)
			{
				tmp_node[j].global = i;
				d_node[i].local = j;
				j++;
			}
		}
		/***************************************************************************************************************/
		h = sd_node;
		k =0;
		/**********************************counting subdomain ghost nodes***********************************************/
		h = sd_node;
		k =0;
		for (i=1; i<=sd_node; i++)
		{
			if (d_node[tmp_node[i].global].corner_ID != 0 || d_node[tmp_node[i].global].loc != 0)
			{
					k++;
			}
			for (j=0; j<=5; j++)
			{
				if (d_node[tmp_node[i].global].n_n[j] != 0 && d_node[tmp_node[i].global].n_n[j] != myrank)
				{
					h = h+6;
				}
			}
		}
		nodes = h;
		sd_gh = h-sd_node;
		
		/***********************************************************************************************/
		tmp_c = (int *)malloc((size+1)*sizeof(int));
		c = (int *)malloc((size+1)*sizeof(int));
		b = (int **)malloc((size+1)*sizeof(int *));
		proc_node = (int *)malloc((size+1)*sizeof(int));
		gb = (int **)malloc((size+1)*sizeof(int *));
		loc_dat = (int **)malloc((size+1)*sizeof(int *));
		recv_b = (int **)malloc((size+1)*sizeof(int *));
		recv_c = (int *)malloc((size+1)*sizeof(int));
		
		for(i=0;i<=size;i++)
		{
			b[i] = (int*)malloc((sd_gh+10)*sizeof(int));
		}
		for(i=0;i<=size;i++)
		{
			gb[i] = (int*)malloc((sd_gh+10)*sizeof(int));
		}
		for(i=0;i<=size;i++)
		{
			loc_dat[i] = (int*)malloc((sd_gh+10)*sizeof(int));
		}
		for(i=0;i<=size;i++)
		{
			recv_b[i] = (int*)malloc((sd_gh+10)*sizeof(int));
		}
		node = (MNODE*)malloc((nodes+10+(k*4))*sizeof(MNODE ));
		singular = (singul*)malloc((nodes+10+(k*4))*sizeof(singul));
		sd_gh_node = (SUB_DOM*)malloc((nodes+10+(k*4))*sizeof(SUB_DOM));
		
		memset(sd_gh_node, 0, sizeof(sd_gh_node));
		
		
		all_bou_node =k;
		temp_boundary_nodes = (int*)malloc((k+10)*sizeof(int));		
		memset(temp_boundary_nodes, 0, sizeof(temp_boundary_nodes));	
		
		for(i=0; i<=nodes+7+(k*4); i++)
		{
			node[i].x=0;
			node[i].y=0;
			node[i].z=0;
			node[i].e[0]=0;
			node[i].e[1]=0;
			node[i].e[2]=0;
			node[i].e[3]=0;	
			node[i].e[4]=0;
			node[i].e[5]=0;
			node[i].e[6]=0;
			node[i].e[7]=0;	
			node[i].ID=0;
			node[i].n_n[0]=0;
			node[i].n_n[1]=0;
			node[i].n_n[2]=0;
			node[i].n_n[3]=0;
			node[i].n_n[4]=0;
			node[i].n_n[5]=0;
		//	node[i].location=0;
		//	node[i].corner=0;
			node[i].corner_ID=0;
			node[i].val=0;
			node[i].loc=0;
		//	node[i].weno=0;
			node[i].proc=0;
			node[i].local=0;
			node[i].global=0;
			node[i].req=0;
			singular[i].n_n[0]=0;
			singular[i].n_n[1]=0;
			singular[i].n_n[2]=0;
			singular[i].n_n[3]=0;
			singular[i].n_n[4]=0;
			singular[i].n_n[5]=0;
		}
		/**********************************************************************************************/
		j = 1;
		for (i=1; i<=NUMNP; i++)
		{
			if (d_node[i].proc == myrank)
			{
				node[j].global = i;
			//	node[j].location = d_node[i].location;
				node[j].loc = d_node[i].loc;
				node[j].val = d_node[i].val;
				node[j].corner_ID = d_node[i].corner_ID;
				d_node[i].local = j;
				j++;
			}
		}
		free(tmp_node);
		tmp_node = NULL;
		/***************GLOBAL AND LOCAL NUMBERING OF NODES AND ITS NEIGHBOURS*************************/
		/**********DOMAIN NEIGHBOUR AND SUBDOMAIN BOUNDARY POINTS PROCS IDENTIFICATION*****************/
		
		h = sd_node+1;
		gh = sd_node;
		k =0;
		temp = 0;
		for (i=1; i<=sd_node; i++)
		{
			//if ((d_node[d_node[node[i].global].n_n[0]].proc == myrank && d_node[d_node[node[i].global].n_n[1]].proc == myrank && d_node[d_node[node[i].global].n_n[2]].proc == myrank && d_node[d_node[node[i].global].n_n[3]].proc == myrank && d_node[d_node[node[i].global].n_n[4]].proc == myrank && d_node[d_node[node[i].global].n_n[5]].proc == myrank) || (d_node[node[i].global].corner_ID != 0 || d_node[node[i].global].loc != 0 ) )
			//{
				node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
				node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
				node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
				node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
				node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
				node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
			//	node[i].location = d_node[node[i].global].location;
				node[i].val = d_node[node[i].global].val;
				node[i].loc = d_node[node[i].global].loc;
				node[i].corner_ID = d_node[node[i].global].corner_ID; 
				node[i].ID = d_node[node[i].global].ID; 
				node[i].proc = d_node[node[i].global].proc; 
				if (d_node[node[i].global].corner_ID != 0 || d_node[node[i].global].loc != 0)
				{
					temp_boundary_nodes[temp] = i;
					temp++;
				}
			//}
		}
		
		h = sd_node+1;
		gh = sd_node;
		k =0;
		for (i=1; i<=sd_node; i++)
		{		
			for (k=0; k <= 5; k++)
			{
				if (d_node[d_node[node[i].global].n_n[k]].proc != myrank && d_node[node[i].global].n_n[k] != 0 && d_node[d_node[node[i].global].n_n[k]].proc > 0 && d_node[d_node[node[i].global].n_n[k]].proc < size)
				{
					if ((k == 0 || k == 1) && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k+2] = i;
						//	node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k+2] = i;
						}
						if (k == 0)
						{
							node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						}
						if (k == 1)
						{
							node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						}
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
						node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}
					
					if ((k == 2 || k == 3) && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k-2] = i;
						//	node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k-2] = i;
						}
						if (k == 2)
						{
							node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						}
						if (k == 3)
						{
							node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						}
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
						node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}
					
					if (k == 4 && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k+1] = i;
						//	node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k+1] = i;
						}
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}
					
					if (k == 5 && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k-1] = i;
					//		node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k-1] = i;
						}
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}					
				}
			
				
				if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[node[i].global].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc < size)
				{
					if ((k == 0 || k == 1) && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k+2] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+2] = d_node[d_node[node[i].global].n_n[k]].local;
						}						
					}
				
					if ((k == 2 || k == 3) && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k-2] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-2] = d_node[d_node[node[i].global].n_n[k]].local;
						}
					}
					
					if (k == 4 && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k+1] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+1] = d_node[d_node[node[i].global].n_n[k]].local;
						}					
					}
				
					if (k == 5 && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k-1] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-1] = d_node[d_node[node[i].global].n_n[k]].local;
						}						
					}
				
				
					if (k == 0 || k == 1)
					{
						if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+2] = h;  
						}
					}
					
					if (k == 2 || k == 3)
					{
						if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-2] = h;  
						}
					}
				
					if (k == 4)
					{
						if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+1] = h;  
						}
					}
				
					if (k == 5)
					{
						if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-1] = h;  
						}
					}
				}

				if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k] != 0  && d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc < size)
				{
					if (k ==0 || k==1)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k+2] = node[node[i].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+2] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						}
					}
					
					if (k ==2 || k==3)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k-2] = node[node[i].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-2] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						}
					}
				
					if (k ==4)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k+1] = node[node[i].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+1] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						}
					}
					
					if (k ==5)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k-1] = node[node[i].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-1] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						}
					}
				
				
					if (k ==0 || k==1)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank )
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+2]= h;
						}
					}	
					
					if (k ==2 || k==3)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-2]= h;
						}
					}
					
					if (k ==4)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank  )
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+1]= h;
						}
					}
				
					if (k ==5)
					{
						if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank  )
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-1]= h;
						}
					}	
				}
				if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc < size)
				{
					if (k == 0 || k == 1)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k+2] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+2] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						}
					}
					
					if (k == 2 || k == 3)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k-2] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-2] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						}
					}
				
					if (k == 4)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k+1] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if(d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+1] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						}
					}
					
					if (k == 5)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
						{
							node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
							node[h].n_n[k-1] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
					//		node[h].location = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
							d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-1] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						}
					}
			
					if (k == 0 || k == 1)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+2]= h;
						}
					}
					
					if (k == 2 || k == 3)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-2]= h;
						}
					}
				
					if (k == 4)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k+1]= h;
						}
					}
					
					if (k == 5)
					{
						if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
						{
							node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k-1]= h;
						}
					}
				}
			}
		}
		
		nodes = h;
		sd_gh = gh-sd_node;
		all_bou_node =temp;
		
		for (i=h-sd_gh; i<=nodes; i++)
		{		
			for (k=0; k <= 5; k++)
			{
				if (d_node[d_node[node[i].global].n_n[k]].proc != myrank && d_node[node[i].global].n_n[k] != 0 && d_node[d_node[node[i].global].n_n[k]].proc > 0 && d_node[d_node[node[i].global].n_n[k]].proc < size)
				{
					if ((k == 0 || k == 1) && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k+2] = i;
						//	node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k+2] = i;
						}
						if (k == 0)
						{
							node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						}
						if (k == 1)
						{
							node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						}
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
						node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}
					
					if ((k == 2 || k == 3) && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k-2] = i;
						//	node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k-2] = i;
						}
						if (k == 2)
						{
							node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						}
						if (k == 3)
						{
							node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						}
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
						node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}
					
					if (k == 4 && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k+1] = i;
						//	node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k+1] = i;
						}
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}
					
					if (k == 5 && node[i].n_n[k] == 0)
					{
						if (d_node[d_node[node[i].global].n_n[k]].local == 0)
						{
							node[i].n_n[k] = h;
							node[h].n_n[k-1] = i;
					//		node[h].location = d_node[d_node[node[i].global].n_n[k]].location;
							node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
							node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
							sd_gh_node[h].global = d_node[node[i].global].n_n[k];
							node[h].global = d_node[node[i].global].n_n[k];
							d_node[d_node[node[i].global].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[node[i].global].n_n[k]].local != 0)
						{
							node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k-1] = i;
						}
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
						node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					//	node[i].location = d_node[node[i].global].location;
						node[i].val = d_node[node[i].global].val;
						node[i].loc = d_node[node[i].global].loc;
						node[i].corner_ID = d_node[node[i].global].corner_ID;
						node[i].ID = d_node[node[i].global].ID;
					}					
				}
			
				
				if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[node[i].global].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc < size)
				{
					if ((k == 0 || k == 1) && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k+2] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+2] = d_node[d_node[node[i].global].n_n[k]].local;
						}						
					}
				
					if ((k == 2 || k == 3) && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k-2] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-2] = d_node[d_node[node[i].global].n_n[k]].local;
						}
					}
					
					if (k == 4 && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k+1] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+1] = d_node[d_node[node[i].global].n_n[k]].local;
						}					
					}
				
					if (k == 5 && node[node[i].n_n[k]].n_n[k] == 0)
					{
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
						{
							node[node[i].n_n[k]].n_n[k] = h;
							node[h].n_n[k-1] = node[i].n_n[k];
					//		node[h].location = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].location;
							node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
							node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
							node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
							node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
							sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
							sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
							d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
							gh++;
							h++;
						}
						if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
						{
							node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
							node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-1] = d_node[d_node[node[i].global].n_n[k]].local;
						}						
					}
				}
				
				if (k == 0 || k == 1)
				{
					if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+2] = d_node[d_node[node[i].global].n_n[k]].local;  
					}
				}
				
				if (k == 2 || k == 3)
				{
					if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-2] = d_node[d_node[node[i].global].n_n[k]].local;  
					}
				}
			
				if (k == 4)
				{
					if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k+1] = d_node[d_node[node[i].global].n_n[k]].local;  
					}
				}
			
				if (k == 5)
				{
					if(d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k-1] = d_node[d_node[node[i].global].n_n[k]].local;  
					}
				}
				
			}
		}
		
		nodes = h;
		sd_gh = gh-sd_node;
		all_bou_node =temp;
		
		
		
		
		
		
		for (i=1; i<=nodes;i++)
		{
			node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
			node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
			node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
			node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
			node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
			node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
			node[i].loc = d_node[node[i].global].loc;
			node[i].x = d_node[node[i].global].x;
			node[i].y = d_node[node[i].global].y;
			node[i].z = d_node[node[i].global].z;
			node[i].val = d_node[node[i].global].val;
		//	node[i].location = d_node[node[i].global].location;
			node[i].corner_ID = d_node[node[i].global].corner_ID;
			node[i].ID = d_node[node[i].global].ID;
 		}
		
		/*****************SEGREGATING NODES OF OTHER PROCESSORS*******************************************/
		for(i=0; i<=size; i++)
		{
			tmp_c[i] = 0;
		}
		for (j=sd_node+1; j< h; j++)
		{
			tmp_c[sd_gh_node[j].proc]++;
		}
		//proc_node = NULL;
		j = 0;
		for (i=1; i<size; i++)
		{
			if(tmp_c[i] != 0 && i != myrank)
			{
				c[j] = i;				
				proc_node[c[j]] = tmp_c[i];
				d[j] = j;
				j++;
			}
		}
		neigh_pro = j;

		/***********************Arranging list of neighbour procs in acending order********************/
		for(i=0; i<neigh_pro; i++)
		{
			for(j=0;j<neigh_pro;j++)
			{
				if(c[j] > c[j+1] && neigh_pro > (j+1))
				{
					temp = c[j+1];
					c[j+1] = c[j];
					c[j] = temp;
					temp = d[j+1];
					d[j+1] = d[j];
					d[j] = temp;
				}
			}
		}
		
/*		for (j=0; j< neigh_pro; j++)
		{	
			k=0;
			for (i=sd_node+1; i<=gh; i++)
			{
				if (sd_gh_node[i].proc == c[j])
				{
					b[j][k] = sd_gh_node[i].global;
					k++;
				}
			}
		}
	*/
		
		MPI_Send(&neigh_pro, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);
		
		position=0;
		MPI_Pack_size(neigh_pro,MPI_INT,MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		for (i=0; i<neigh_pro; i++)
		{
			MPI_Pack(&c[i],1,MPI_INT,buffer,memsize,&position, MPI_COMM_WORLD);
		}
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		free(buffer);
		buffer = NULL;
		

		
		int *all_neigh_pro, **all_c_list;
		all_neigh_pro = (int *)malloc((size+1)*sizeof(int));
		all_c_list = (int **)malloc((size+1)*sizeof(int *));
		
		for (i=1; i<=size-1; i++)
		{			
			MPI_Recv(&all_neigh_pro[i], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);		
			
			all_c_list[i] = (int *)malloc((all_neigh_pro[i]+1)*sizeof(int));
			position = 0;
			MPI_Pack_size(all_neigh_pro[i], MPI_INT,MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);
			MPI_Recv(buffer, memsize, MPI_PACKED, 0, myrank, MPI_COMM_WORLD, &status);		
			for (j=1; j<=all_neigh_pro[i]; j++)
			{
				MPI_Unpack(buffer,memsize,&position,&glob,1,MPI_INT,MPI_COMM_WORLD);
				all_c_list[i][j] = glob;
			}
			free(buffer);
			buffer = NULL;
		}
		
		for (i=1;i<=size-1;i++)
		{
			for (j=1;j<=all_neigh_pro[i];j++)
			{
				if (all_c_list[i][j] == myrank && proc_node[i] == 0 && i != myrank)
				{
					neigh_pro++;
					c[neigh_pro-1] = i;
					proc_node[i] = 0;
				}
			}
		}
		
		
		for(i=0; i<neigh_pro; i++)
		{
			for(j=0;j<neigh_pro;j++)
			{
				if(c[j] > c[j+1] && neigh_pro > (j+1))
				{
					temp = c[j+1];
					c[j+1] = c[j];
					c[j] = temp;
					temp = d[j+1];
					d[j+1] = d[j];
					d[j] = temp;
				}
			}
		}		
		
		for (j=0; j< neigh_pro; j++)
		{	
			k=0;
			for (i=sd_node+1; i<=gh; i++)
			{
				if (sd_gh_node[i].proc == c[j])
				{
					b[j][k] = sd_gh_node[i].global;
					k++;
				}
			}
		}
		
		/****************************data send to master for proc request checking*********************/
		temp =0;
		for (i=0; i<neigh_pro; i++)
		{
			temp = temp+proc_node[c[i]];
		}
		MPI_Send(&sd_gh, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);
		position=0;
		MPI_Pack_size(temp,MPI_INT,MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		for (i=0; i<neigh_pro; i++)
		{
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Pack(&b[i][j],1,MPI_INT,buffer,memsize,&position, MPI_COMM_WORLD);
			}
		}
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		free(buffer);
		buffer = NULL;

		/**********************SEND RECIEVE LOCAL NODE NEIGHBOUR DATA**********************************/		
		position=0;
		no_of_tip_recv = 0;
		k = 0;
		for (i=0; i<neigh_pro; i++)
		{
			position=0;
			MPI_Sendrecv(&proc_node[c[i]], 1, MPI_INT, c[i], c[i], &recv_c[c[i]],1, MPI_INT, c[i], myrank, MPI_COMM_WORLD, &status);
			MPI_Pack_size(proc_node[c[i]],MPI_INT,MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);                                            
			MPI_Pack_size(recv_c[c[i]],MPI_INT,MPI_COMM_WORLD,&memsize1);
			buffer2 = malloc(memsize1);                              /***********carefull with buffer1 and buffer2******************/
			position = 0;
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Pack(&b[i][j],1,MPI_INT,buffer,memsize,&position, MPI_COMM_WORLD);
				recv_b[i][j] = d_node[b[i][j]].local;
				if (d_node[b[i][j]].loc >= 27 && d_node[b[i][j]].loc <= 32)
				{
					tip_recv[++k].node = d_node[b[i][j]].local;
					no_of_tip_recv++;
				}
			}
			MPI_Sendrecv(buffer,position,MPI_PACKED,c[i],c[i],buffer2,memsize1,MPI_PACKED,c[i],myrank, MPI_COMM_WORLD, &status);
			position = 0;
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Unpack(buffer2,memsize1,&position,&glob,1,MPI_INT,MPI_COMM_WORLD);
				gb[i][j] = glob;    /*************Neigh processor's needs local node number of nodes with these global numbers frm current procs ***********/
			}
			free(buffer);
			buffer=NULL;
			free(buffer2);
			buffer2=NULL;
		}
		position=0;
		no_of_tip_send = 0;
		k = 0;
		for (i=0;i<neigh_pro;i++)
		{
			for(j=0; j<recv_c[c[i]]; j++)
			{
				loc_dat[i][j] = d_node[gb[i][j]].local;    /************List of node numbers to be sent to neighbour processor********/            
				if (d_node[gb[i][j]].loc >= 27 && d_node[gb[i][j]].loc <= 32)
				{
					no_of_tip_send++;
					tip[++k].i = i;
					tip[k].j = j;
					tip[k].node = d_node[gb[i][j]].local;
					tip[k].proc = d_node[gb[i][j]].proc;
				}  
			}
		}
		
		for (i=0; i<neigh_pro; i++)
		{
			MPI_Pack_size(recv_c[c[i]],MPI_INT,MPI_COMM_WORLD,&memsize1);
			buffer2 = malloc(memsize1);                                          /***********carefull with buffer1 and buffer2******************/
			MPI_Pack_size(proc_node[c[i]],MPI_INT,MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);
			position = 0;
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&loc_dat[i][j],1,MPI_INT,buffer2,memsize1,&position, MPI_COMM_WORLD);
			}
			
			MPI_Sendrecv(buffer2,memsize1,MPI_PACKED,c[i],c[i],buffer,memsize,MPI_PACKED,c[i],myrank, MPI_COMM_WORLD, &status);
			position = 0;
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&loca,1,MPI_INT,MPI_COMM_WORLD);
				b[i][j] = loca;
			}
			free(buffer);
			buffer=NULL;
			free(buffer2);
			buffer2=NULL;
		}
		
		
		/**************************Obtaining boundary nodes data from master to slave*************************/
		sprintf(filename,"inlet_%d.txt",myrank);
		fnode=fopen(filename,"rt");
		inlet_node = (int*)malloc((inl_node+1)*sizeof(int));
		i=0;
		while(fgets(line, 200, fnode) != NULL)
		{
			sscanf(line,"%d",&temp);
			inlet_node[i] = d_node[temp].local;
			i++;
		}	
		//temp_inl_node = i;
		fclose(fnode);
		
		sprintf(filename,"outlet_%d.txt",myrank);
		fnode=fopen(filename,"rt");
		outlet_node = (int*)malloc((out_node+1)*sizeof(int));
		i=0;
		while(fgets(line, 200, fnode) != NULL)
		{
			sscanf(line,"%d",&temp);
			outlet_node[i] = d_node[temp].local;
			i++;
		}	
		//temp_out_node = i;
		fclose(fnode);
		
		sprintf(filename,"boundary_%d.txt",myrank);
		fnode=fopen(filename,"rt");
		boundary_node = (int*)malloc((bou_node+1)*sizeof(int));
		i=0;
		while(fgets(line, 200, fnode) != NULL)
		{
			sscanf(line,"%d",&temp);
			boundary_node[i] = d_node[temp].local;
			i++;
		}
		//temp_bou_node = i;
		fclose(fnode);
		
		sprintf(filename,"wall_%d.txt",myrank);
		fnode=fopen(filename,"rt");
		wall_node = (int*)malloc((wal_node+1)*sizeof(int));
		i=0;
		while(fgets(line, 200, fnode) != NULL)
		{
			sscanf(line,"%d",&temp);
			wall_node[i] = d_node[temp].local;
			i++;
		}
		//temp_wal_node = i;
		fclose(fnode);
		/************************sub_domain boundary conditions*****************************/
		sd_inl_node =0;
		sd_out_node =0;
		sd_wal_node =0;
		sd_bou_node =0;
		for (i = sd_node+1; i<nodes;i++)
		{
			if(node[i].ID == 1)
			{
				sd_inl_node++;
			}
			if(node[i].ID == 2)
			{
				sd_out_node++;
			}
			if(node[i].ID == 3 || node[i].ID == 10)
			{
				sd_wal_node++;
			}
			if(node[i].ID == 4)
			{
				sd_bou_node++;
			}
		}
		
		sd_inlet_node = (int *)malloc((sd_inl_node+1)*sizeof(int));
		sd_outlet_node = (int *)malloc((sd_out_node+1)*sizeof(int));
		sd_wall_node = (int *)malloc((sd_wal_node+1)*sizeof(int));
		sd_boundary_node = (int *)malloc((sd_bou_node+1)*sizeof(int));
		temp = 0;
		for (i = sd_node+1; i<=nodes;i++)
		{
			if(node[i].ID == 1)
			{
				sd_inlet_node[temp++] = i;
			}			
		}
		temp = 0;
		for (i = sd_node+1; i<=nodes;i++)
		{
			if(node[i].ID == 2)
			{
				sd_outlet_node[temp++] = i;
			}			
		}
		temp = 0;
		for (i = sd_node+1; i<=nodes;i++)
		{
			if(node[i].ID == 3 || node[i].ID == 10)
			{
				sd_wall_node[temp++] = i;
			}			
		}
		temp = 0;
		for (i = sd_node+1; i<=nodes;i++)
		{
			if(node[i].ID == 4)
			{
				sd_boundary_node[temp++] = i;
			}			
		}
		/*************************************************************************************/
/*		int *temp_all_boundary_nodes, *temp_saving, temp_bound_count;
		temp_all_boundary_nodes = (int *)malloc((sd_node+sd_inl_node+sd_out_node+sd_wal_node+sd_bou_node)*sizeof(int));
		temp_saving = (int *)malloc((sd_node+sd_inl_node+sd_out_node+sd_wal_node+sd_bou_node)*sizeof(int));
		
		temp = all_bou_node;
		all_bou_node = all_bou_node+sd_inl_node+sd_out_node+sd_wal_node+sd_bou_node;
		all_boundary_nodes = (int*)malloc((all_bou_node+k+10)*sizeof(int));		
		memset(all_boundary_nodes, 0, sizeof(all_boundary_nodes));
		temp_bound_count =0;
		for (i=0; i<sd_node+sd_inl_node+sd_out_node+sd_wal_node+sd_bou_node; i++)
		{
			temp_saving[i] = 0;
			temp_all_boundary_nodes[i] = 0;
		}
		for (i=0; i<temp; i++)
		{
			if (temp_saving[temp_boundary_nodes[i]] == 0)
			{
				all_boundary_nodes[i] = temp_boundary_nodes[i];
				temp_saving[temp_boundary_nodes[i]] = 1;
				temp_bound_count++;
			}
		}	
		//temp++;
		for (i=0; i<sd_inl_node; i++)
		{
			if (temp_saving[sd_inlet_node[i]] == 0)
			{
				all_boundary_nodes[temp++] = sd_inlet_node[i];
				temp_saving[sd_inlet_node[i]] = 1;
				temp_bound_count++;
			}
		}
		for (i=0; i<sd_out_node; i++)
		{
			if (temp_saving[sd_outlet_node[i]] == 0)
			{
				all_boundary_nodes[temp++] = sd_outlet_node[i];
				temp_saving[sd_outlet_node[i]] = 1;
				temp_bound_count++;
			}
		}
		for (i=0; i<sd_wal_node; i++)
		{
			if (temp_saving[sd_wall_node[i]] == 0)
			{
				all_boundary_nodes[temp++] = sd_wall_node[i];
				temp_saving[sd_wall_node[i]] = 1;
				temp_bound_count++;
			}
		}
		for (i=0; i<sd_bou_node; i++)
		{
			if (temp_saving[sd_boundary_node[i]] == 0)
			{
				all_boundary_nodes[temp++] = sd_boundary_node[i];
				temp_saving[sd_boundary_node[i]] = 1;
				temp_bound_count++;
			}
		}
		all_bou_node = temp_bound_count;
		
		free(temp_all_boundary_nodes);
		temp_all_boundary_nodes = NULL;
		free(temp_saving);
		temp_saving = NULL;
		
		*/
		
		temp = all_bou_node;
		all_bou_node = all_bou_node+sd_inl_node+sd_out_node+sd_wal_node+sd_bou_node;
		all_boundary_nodes = (int*)malloc((all_bou_node+k+10)*sizeof(int));		
		memset(all_boundary_nodes, 0, sizeof(all_boundary_nodes));
		for (i=0; i<all_bou_node; i++)
		{
			all_boundary_nodes[i] = temp_boundary_nodes[i];
		}		
		for (i=0; i<sd_inl_node; i++)
		{
			all_boundary_nodes[temp++] = sd_inlet_node[i];
		}
		for (i=0; i<sd_out_node; i++)
		{
			all_boundary_nodes[temp++] = sd_outlet_node[i];
		}
		for (i=0; i<sd_wal_node; i++)
		{
			all_boundary_nodes[temp++] = sd_wall_node[i];
		}
		for (i=0; i<sd_bou_node; i++)
		{
			all_boundary_nodes[temp++] = sd_boundary_node[i];
		}
		
		/*************************************************************************************/
		metric = (TRANSFORMED *)malloc((nodes+10+(all_bou_node*5))*sizeof(TRANSFORMED));
		jacobian = (JACOB *)malloc((nodes+10+(all_bou_node*5))*sizeof(JACOB));
		deter = (DETERM *)malloc((nodes+10+(all_bou_node*5))*sizeof(DETERM));
		det = (double *)malloc((nodes+10+(all_bou_node*5))*sizeof(double));
		Q = (Qip **)malloc((nodes+10+(all_bou_node*5))*sizeof(Qip*));
		for (i=0; i<(nodes+10+(all_bou_node*5)); i++)
		{
			Q[i] = (Qip *)malloc(10*sizeof(Qip));
		}
		
		t_metric = (TRANSFORMED *)malloc((NUMNP+10)*sizeof(TRANSFORMED));
		t_jacobian = (JACOB *)malloc((NUMNP+10)*sizeof(JACOB));
		t_det = (double *)malloc((NUMNP+10)*sizeof(double));
		
		memset(jacobian, 0, sizeof(jacobian));
		memset(deter, 0, sizeof(deter));
		memset(metric, 0, sizeof(metric));
		memset(det,0, sizeof(det));
		memset(t_jacobian, 0, sizeof(t_jacobian));
		memset(t_metric, 0, sizeof(t_metric));
		memset(t_det,0, sizeof(t_det));
		
			
		MPI_Pack_size((NUMNP+10)*5, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].x_zeta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].x_eta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].x_xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].y_zeta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].y_eta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		}
		free(Bcast_buffer3);
		Bcast_buffer3 = NULL;	
		
		MPI_Pack_size((NUMNP+10)*4, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].y_xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].z_zeta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].z_eta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].z_xi, 1, MPI_DOUBLE, MPI_COMM_WORLD); 
		}
		free(Bcast_buffer4);
		Bcast_buffer4 = NULL;	
		
		
		for(i=1;i<=NUMNP;i++)
		{
			t_det[i] = (((t_jacobian[i].x_zeta)*((t_jacobian[i].y_eta*t_jacobian[i].z_xi)-(t_jacobian[i].y_xi*t_jacobian[i].z_eta)))+((t_jacobian[i].x_eta)*((t_jacobian[i].y_xi*t_jacobian[i].z_zeta)-(t_jacobian[i].y_zeta*t_jacobian[i].z_xi)))+((t_jacobian[i].x_xi)*((t_jacobian[i].y_zeta*t_jacobian[i].z_eta)-(t_jacobian[i].y_eta*t_jacobian[i].z_zeta))));
			t_metric[i].zeta_x = (1.0/t_det[i])*(t_jacobian[i].y_eta*t_jacobian[i].z_xi-t_jacobian[i].y_xi*t_jacobian[i].z_eta);
			t_metric[i].zeta_y = (1.0/t_det[i])*(t_jacobian[i].z_eta*t_jacobian[i].x_xi-t_jacobian[i].z_xi*t_jacobian[i].x_eta);
			t_metric[i].zeta_z = (1.0/t_det[i])*(t_jacobian[i].x_eta*t_jacobian[i].y_xi-t_jacobian[i].x_xi*t_jacobian[i].y_eta);
			t_metric[i].eta_x = (1.0/t_det[i])*(t_jacobian[i].y_xi*t_jacobian[i].z_zeta-t_jacobian[i].y_zeta*t_jacobian[i].z_xi);
			t_metric[i].eta_y = (1.0/t_det[i])*(t_jacobian[i].z_xi*t_jacobian[i].x_zeta-t_jacobian[i].z_zeta*t_jacobian[i].x_xi);
			t_metric[i].eta_z = (1.0/t_det[i])*(t_jacobian[i].x_xi*t_jacobian[i].y_zeta-t_jacobian[i].x_zeta*t_jacobian[i].y_xi);		
			t_metric[i].xi_x = (1.0/t_det[i])*(t_jacobian[i].y_zeta*t_jacobian[i].z_eta-t_jacobian[i].y_eta*t_jacobian[i].z_zeta);
			t_metric[i].xi_y = (1.0/t_det[i])*(t_jacobian[i].z_zeta*t_jacobian[i].x_eta-t_jacobian[i].z_eta*t_jacobian[i].x_zeta);
			t_metric[i].xi_z = (1.0/t_det[i])*(t_jacobian[i].x_zeta*t_jacobian[i].y_eta-t_jacobian[i].x_eta*t_jacobian[i].y_zeta);	
		
		}
		
		/****************************READING TRANSFORMATION DATA***********************************/
		i =1;
/*		fnode = fopen("transformation_data.txt","rt");
		if (fnode != NULL)
		{
			while(fgets(line, 400, fnode) != NULL)
			{
				sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t_jacobian[i].x_zeta, &t_jacobian[i].x_eta, &t_jacobian[i].x_xi, &t_jacobian[i].y_zeta, &t_jacobian[i].y_eta, &t_jacobian[i].y_xi, &t_jacobian[i].z_zeta, &t_jacobian[i].z_eta, &t_jacobian[i].z_xi, &t_metric[i].zeta_x, &t_metric[i].eta_x, &t_metric[i].xi_x, &t_metric[i].zeta_y, &t_metric[i].eta_y, &t_metric[i].xi_y, &t_metric[i].zeta_z, &t_metric[i].eta_z, &t_metric[i].xi_z, &t_det[i]);
				i++;
			}	
			fclose(fnode);
		}
		*/
		for (i=1; i<=nodes; i++)
		{
			jacobian[i].x_zeta = t_jacobian[node[i].global].x_zeta;
			jacobian[i].x_eta = t_jacobian[node[i].global].x_eta;
			jacobian[i].x_xi = t_jacobian[node[i].global].x_xi;
			jacobian[i].y_zeta = t_jacobian[node[i].global].y_zeta;
			jacobian[i].y_eta = t_jacobian[node[i].global].y_eta;
			jacobian[i].y_xi = t_jacobian[node[i].global].y_xi;
			jacobian[i].z_zeta = t_jacobian[node[i].global].z_zeta;
			jacobian[i].z_eta = t_jacobian[node[i].global].z_eta;
			jacobian[i].z_xi = t_jacobian[node[i].global].z_xi;
			metric[i].zeta_x = t_metric[node[i].global].zeta_x;
			metric[i].eta_x = t_metric[node[i].global].eta_x;
			metric[i].xi_x = t_metric[node[i].global].xi_x;
			metric[i].zeta_y = t_metric[node[i].global].zeta_y;
			metric[i].eta_y = t_metric[node[i].global].eta_y;
			metric[i].xi_y = t_metric[node[i].global].xi_y;
			metric[i].zeta_z = t_metric[node[i].global].zeta_z;
			metric[i].eta_z = t_metric[node[i].global].eta_z;
			metric[i].xi_z = t_metric[node[i].global].xi_z;
			det[i] = t_det[node[i].global];
		}
		/****************************Free arrays******************************************************/
		free(t_jacobian);
		t_jacobian = NULL;
		free(t_metric);
		t_metric = NULL;
		free(t_det);
		t_det = NULL;
		free(temp_boundary_nodes);
		temp_boundary_nodes = NULL;
		free(sd_gh_node);
		sd_gh_node = NULL;
		/**************************Memory Reallocation for ghost points*******************************/
		g_node=sd_node+sd_gh+1;		
		/*********************ghost elements and nodes creation***************************************/

		for(i=0; i<all_bou_node; i++)
		{			
			if (node[all_boundary_nodes[i]].corner_ID == 0)
			{
				if (node[all_boundary_nodes[i]].loc == 4)
				{
					j = 3; /********no node on WEST**************/
					k = 1;
					m = 4;	
					gar5 = 1;						
				}
			
				if (node[all_boundary_nodes[i]].loc == 3)
				{
					j = 2; /********no element on SOUTH**************/
					k = 0;
					m = 3;
					gar5 = 2;			
				}
				
				if (node[all_boundary_nodes[i]].loc == 2)
				{
					j = 1; /********no element on EAST**************/
					k = 3;
					m = 2;
					gar5 = 1;
				}
				
				if (node[all_boundary_nodes[i]].loc == 1)
				{
					j = 0; /********no element on NORTH**************/
					k = 2;
					m = 1;		
					gar5 = 2;	
				}
				
				if (node[all_boundary_nodes[i]].loc == 6)
				{
					j = 4; /********no element on NORTH**************/
					k = 5;
					m = 1;		
					gar5 = 2;	
				}
				
				if (node[all_boundary_nodes[i]].loc == 5)
				{
					j = 5; /********no element on NORTH**************/
					k = 4;
					m = 1;		
					gar5 = 2;	
				}
				
				if (node[all_boundary_nodes[i]].loc > 14 && node[all_boundary_nodes[i]].loc < 27)
				{
					gar5 = 4;
				}
				
				if (node[all_boundary_nodes[i]].loc >= 27 && node[all_boundary_nodes[i]].loc <= 32)
				{
					gar5 = 3;
				}
				
				switch (gar5)
				{
					case 1:
						node[all_boundary_nodes[i]].n_n[j] = g_node;
						node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
						det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
					
						node[g_node-1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node-1;
						metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
						
						node[g_node-1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node-1;
						metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
					break;
			
					case 2:
						node[all_boundary_nodes[i]].n_n[j] = g_node;
						node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
						det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
					
						node[g_node-1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node-1;
						metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
						
						node[g_node-1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node-1;
						metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
					break;
					
					case 3:
						if (node[all_boundary_nodes[i]].loc == 27 || node[all_boundary_nodes[i]].loc == 28)
						{
							temp6 = 3;
						}
						if (node[all_boundary_nodes[i]].loc == 29)
						{
							temp6 = 4;
						}	
						if (node[all_boundary_nodes[i]].loc == 30)
						{
							temp6 = 4;
						}
						
						for (temp = 0; temp<temp6; temp++)
						{
							if (temp == 0 && (node[all_boundary_nodes[i]].loc == 28 || node[all_boundary_nodes[i]].loc == 29 || node[all_boundary_nodes[i]].loc == 30))
							{
								j = 0; /********no element on NORTH**************/
								k = 2;
								m = 1;	
							}
							if (temp == 1 && (node[all_boundary_nodes[i]].loc == 28 || node[all_boundary_nodes[i]].loc == 29 || node[all_boundary_nodes[i]].loc == 30))
							{
								j = 2; /********no element on SOUTH**************/
								k = 0;
								m = 3;
							}
							if (temp == 0 && node[all_boundary_nodes[i]].loc == 27)
							{
								j = 4; /********no element on NORTH**************/
								k = 5;
								m = 1;	
							}
							if (temp == 1 && node[all_boundary_nodes[i]].loc == 27)
							{
								j = 5; /********no element on SOUTH**************/
								k = 4;
								m = 3;
							}
							if (temp == 2)
							{
								j = 1; /********no element on EAST**************/
								k = 3;
								m = 2;
							}
							if (temp == 3 && node[all_boundary_nodes[i]].loc == 29)
							{
								j = 5;
								k = 4;
							}
							if (temp == 3 && node[all_boundary_nodes[i]].loc == 30)
							{
								j = 4;
								k = 5;
							}
							if (temp <=1 || temp == 3)
							{
								if (temp == 0)
								{
									singular[all_boundary_nodes[i]].n_n[j] = node[all_boundary_nodes[i]].n_n[j];
									singular[all_boundary_nodes[i]].n_n[k] = node[all_boundary_nodes[i]].n_n[k];
								}
								node[all_boundary_nodes[i]].n_n[j] = g_node;
								node[g_node].n_n[k] = all_boundary_nodes[i];
								metric[g_node].zeta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_x;
								metric[g_node].eta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_x;
								metric[g_node].xi_x = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_x;
								metric[g_node].zeta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_y;
								metric[g_node].eta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_y;
								metric[g_node].xi_y = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_y;
								metric[g_node].zeta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_z;
								metric[g_node].eta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_z;
								metric[g_node].xi_z = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_z;
								det[g_node] = det[singular[all_boundary_nodes[i]].n_n[k]];
								node[g_node].loc = 100;
								g_node++;
								g_elem++;
						
								node[g_node-1].n_n[j] = g_node;
								node[g_node].n_n[k] = g_node-1;
								metric[g_node].zeta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
								metric[g_node].eta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
								metric[g_node].xi_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
								metric[g_node].zeta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
								metric[g_node].eta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
								metric[g_node].xi_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
								metric[g_node].zeta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
								metric[g_node].eta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
								metric[g_node].xi_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
								det[g_node] = det[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]];
								node[g_node].loc = 100;
								g_node++;
								g_elem++;
						
								node[g_node-1].n_n[j] = g_node;
								node[g_node].n_n[k] = g_node-1;
								metric[g_node].zeta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
								metric[g_node].eta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
								metric[g_node].xi_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
								metric[g_node].zeta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
								metric[g_node].eta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
								metric[g_node].xi_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
								metric[g_node].zeta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
								metric[g_node].eta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
								metric[g_node].xi_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
								det[g_node] = det[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
								node[g_node].loc = 100;
								g_node++;
								g_elem++;
							}
							if (temp == 2)
							{
								singular[all_boundary_nodes[i]].n_n[j] = node[all_boundary_nodes[i]].n_n[j];
								singular[all_boundary_nodes[i]].n_n[k] = node[all_boundary_nodes[i]].n_n[k];
								node[all_boundary_nodes[i]].n_n[j] = g_node;
								node[g_node].n_n[k] = all_boundary_nodes[i];
								metric[g_node].zeta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_x;
								metric[g_node].eta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_x;
								metric[g_node].xi_x = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_x;
								metric[g_node].zeta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_y;
								metric[g_node].eta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_y;
								metric[g_node].xi_y = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_y;
								metric[g_node].zeta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_z;
								metric[g_node].eta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_z;
								metric[g_node].xi_z = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_z;
								det[g_node] = det[singular[all_boundary_nodes[i]].n_n[k]];
								node[g_node].loc = 100;
								g_node++;
								g_elem++;
						
								node[g_node-1].n_n[j] = g_node;
								node[g_node].n_n[k] = g_node-1;
								metric[g_node].zeta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
								metric[g_node].eta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
								metric[g_node].xi_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
								metric[g_node].zeta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
								metric[g_node].eta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
								metric[g_node].xi_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
								metric[g_node].zeta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
								metric[g_node].eta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
								metric[g_node].xi_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
								det[g_node] = det[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]];
								node[g_node].loc = 100;
								g_node++;
								g_elem++;
						
								node[g_node-1].n_n[j] = g_node;
								node[g_node].n_n[k] = g_node-1;
								metric[g_node].zeta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
								metric[g_node].eta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
								metric[g_node].xi_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
								metric[g_node].zeta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
								metric[g_node].eta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
								metric[g_node].xi_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
								metric[g_node].zeta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
								metric[g_node].eta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
								metric[g_node].xi_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
								det[g_node] = det[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
								node[g_node].loc = 100;
								g_node++;
								g_elem++;
							}
						}						
					break;
				
					case 4:
						temp1 =0;
						temp2 =0;
						temp3 =0;
						temp4 =0;
						temp5 =0;
						temp6 =0;
						for (temp = 0; temp<2; temp++)
						{
							temp7 = 0;
							if (node[all_boundary_nodes[i]].n_n[3] == 0 && temp1 == 0 && temp7 == 0)
							{
								j = 3; /********no node on WEST**************/
								k = 1;
								m = 4;						
								temp1++;
								temp7++;
							}
							else if (node[all_boundary_nodes[i]].n_n[2] == 0 && temp2 == 0 && temp7 == 0)
							{
								j = 2; /********no element on SOUTH**************/
								k = 0;
								m = 3;	
								temp2++;
								temp7++;
							}
							else if (node[all_boundary_nodes[i]].n_n[1] == 0 && temp3 == 0 && temp7 == 0)
							{
								j = 1; /********no element on EAST**************/
								k = 3;
								m = 2;
								temp3++;
								temp7++;
							}
							else if (node[all_boundary_nodes[i]].n_n[0] == 0 && temp4 == 0 && temp7 == 0)
							{
								j = 0; /********no element on NORTH**************/
								k = 2;
								m = 1;		
								temp4++;
								temp7++;
							}
							else if (node[all_boundary_nodes[i]].n_n[4] == 0 && temp5 == 0 && temp7 == 0)
							{
								j = 4; /********no element on Z-plus**************/
								k = 5;
								m = 1;		
								temp5++;
								temp7++;
							}
							else if (node[all_boundary_nodes[i]].n_n[5] == 0 && temp6 == 0 && temp7 == 0)
							{
								j = 5; /********no element on Z-minus**************/
								k = 4;
								m = 1;		
								temp6++;
								temp7++;
							}
							
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						
							node[g_node-1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node-1;
							metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
							
							node[g_node-1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node-1;
							metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						}
					break;
				}										
			}
			
			else if (node[all_boundary_nodes[i]].corner_ID != 0)
			{
				for (gar=0; gar<3; gar++)
				{
					if (node[all_boundary_nodes[i]].corner_ID == 1)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 2)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 3)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 4)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;	
						}
					}
					
					if (node[all_boundary_nodes[i]].corner_ID == 5)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 6)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 7)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 8)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;	
						}
					}
					
					switch (gar5)
					{
						case 1:
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						
							node[g_node-1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node-1;
							metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
							
							node[g_node-1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node-1;
							metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						break;

						case 2:
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						
							node[g_node-1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node-1;
							metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
							
							node[g_node-1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node-1;
							metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						break;
					}	
				}
			}
		}

		nodes = g_node;
		/**************************************evaluating metrics at half points***********************************/	
		for(i=1; i<=sd_node; i++)
		{
			metric[i].zeta_xip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_x)-25.0*(metric[node[i].n_n[3]].zeta_x)+150.0*(metric[i].zeta_x)+150.0*(metric[node[i].n_n[1]].zeta_x)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_x)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].zeta_x));
			metric[i].zeta_yip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_y)-25.0*(metric[node[i].n_n[3]].zeta_y)+150.0*(metric[i].zeta_y)+150.0*(metric[node[i].n_n[1]].zeta_y)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_y)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].zeta_y));
			metric[i].zeta_zip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_z)-25.0*(metric[node[i].n_n[3]].zeta_z)+150.0*(metric[i].zeta_z)+150.0*(metric[node[i].n_n[1]].zeta_z)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_z)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].zeta_z));
			metric[i].eta_xip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_x)-25.0*(metric[node[i].n_n[3]].eta_x)+150.0*(metric[i].eta_x)+150.0*(metric[node[i].n_n[1]].eta_x)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_x)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_x));
			metric[i].eta_yip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_y)-25.0*(metric[node[i].n_n[3]].eta_y)+150.0*(metric[i].eta_y)+150.0*(metric[node[i].n_n[1]].eta_y)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_y)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_y));
			metric[i].eta_zip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_z)-25.0*(metric[node[i].n_n[3]].eta_z)+150.0*(metric[i].eta_z)+150.0*(metric[node[i].n_n[1]].eta_z)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_z)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_z));
			metric[i].xi_xip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_x)-25.0*(metric[node[i].n_n[3]].xi_x)+150.0*(metric[i].xi_x)+150.0*(metric[node[i].n_n[1]].xi_x)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_x)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_x));
			metric[i].xi_yip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_y)-25.0*(metric[node[i].n_n[3]].xi_y)+150.0*(metric[i].xi_y)+150.0*(metric[node[i].n_n[1]].xi_y)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_y)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_y));
			metric[i].xi_zip = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_z)-25.0*(metric[node[i].n_n[3]].xi_z)+150.0*(metric[i].xi_z)+150.0*(metric[node[i].n_n[1]].xi_z)-25.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_z)+3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_z));
			deter[i].ip = (1.0/256.0)*(3.0*(det[node[node[i].n_n[3]].n_n[3]])-25.0*(det[node[i].n_n[3]])+150.0*(det[i])+150.0*(det[node[i].n_n[1]])-25.0*(det[node[node[i].n_n[1]].n_n[1]])+3.0*(det[node[node[node[i].n_n[1]].n_n[1]].n_n[1]]));
			
			metric[i].zeta_xim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].zeta_x)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_x)+150.0*(metric[node[i].n_n[3]].zeta_x)+150.0*(metric[i].zeta_x)-25.0*(metric[node[i].n_n[1]].zeta_x)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_x));
			metric[i].zeta_yim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].zeta_y)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_y)+150.0*(metric[node[i].n_n[3]].zeta_y)+150.0*(metric[i].zeta_y)-25.0*(metric[node[i].n_n[1]].zeta_y)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_y));
			metric[i].zeta_zim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].zeta_z)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_z)+150.0*(metric[node[i].n_n[3]].zeta_z)+150.0*(metric[i].zeta_z)-25.0*(metric[node[i].n_n[1]].zeta_z)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_z));
			metric[i].eta_xim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_x)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_x)+150.0*(metric[node[i].n_n[3]].eta_x)+150.0*(metric[i].eta_x)-25.0*(metric[node[i].n_n[1]].eta_x)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_x));
			metric[i].eta_yim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_y)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_y)+150.0*(metric[node[i].n_n[3]].eta_y)+150.0*(metric[i].eta_y)-25.0*(metric[node[i].n_n[1]].eta_y)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_y));	
			metric[i].eta_zim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_z)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_z)+150.0*(metric[node[i].n_n[3]].eta_z)+150.0*(metric[i].eta_z)-25.0*(metric[node[i].n_n[1]].eta_z)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_z));	
			metric[i].xi_xim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_x)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_x)+150.0*(metric[node[i].n_n[3]].xi_x)+150.0*(metric[i].xi_x)-25.0*(metric[node[i].n_n[1]].xi_x)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_x));
			metric[i].xi_yim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_y)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_y)+150.0*(metric[node[i].n_n[3]].xi_y)+150.0*(metric[i].xi_y)-25.0*(metric[node[i].n_n[1]].xi_y)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_y));	
			metric[i].xi_zim = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_z)-25.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_z)+150.0*(metric[node[i].n_n[3]].xi_z)+150.0*(metric[i].xi_z)-25.0*(metric[node[i].n_n[1]].xi_z)+3.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_z));	
			deter[i].im = (1.0/256.0)*(3.0*(det[node[node[node[i].n_n[3]].n_n[3]].n_n[3]])-25.0*(det[node[node[i].n_n[3]].n_n[3]])+150.0*(det[node[i].n_n[3]])+150.0*(det[i])-25.0*(det[node[i].n_n[1]])+3.0*(det[node[node[i].n_n[1]].n_n[1]]));	

			metric[i].zeta_xjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_x)-25.0*(metric[node[i].n_n[2]].zeta_x)+150.0*(metric[i].zeta_x)+150.0*(metric[node[i].n_n[0]].zeta_x)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_x)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_x));
			metric[i].zeta_yjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_y)-25.0*(metric[node[i].n_n[2]].zeta_y)+150.0*(metric[i].zeta_y)+150.0*(metric[node[i].n_n[0]].zeta_y)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_y)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_y));
			metric[i].zeta_zjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_z)-25.0*(metric[node[i].n_n[2]].zeta_z)+150.0*(metric[i].zeta_z)+150.0*(metric[node[i].n_n[0]].zeta_z)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_z)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_z));
			metric[i].eta_xjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_x)-25.0*(metric[node[i].n_n[2]].eta_x)+150.0*(metric[i].eta_x)+150.0*(metric[node[i].n_n[0]].eta_x)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_x)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].eta_x));
			metric[i].eta_yjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_y)-25.0*(metric[node[i].n_n[2]].eta_y)+150.0*(metric[i].eta_y)+150.0*(metric[node[i].n_n[0]].eta_y)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_y)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].eta_y));
			metric[i].eta_zjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_z)-25.0*(metric[node[i].n_n[2]].eta_z)+150.0*(metric[i].eta_z)+150.0*(metric[node[i].n_n[0]].eta_z)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_z)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].eta_z));
			metric[i].xi_xjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_x)-25.0*(metric[node[i].n_n[2]].xi_x)+150.0*(metric[i].xi_x)+150.0*(metric[node[i].n_n[0]].xi_x)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_x)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_x));
			metric[i].xi_yjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_y)-25.0*(metric[node[i].n_n[2]].xi_y)+150.0*(metric[i].xi_y)+150.0*(metric[node[i].n_n[0]].xi_y)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_y)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_y));
			metric[i].xi_zjp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_z)-25.0*(metric[node[i].n_n[2]].xi_z)+150.0*(metric[i].xi_z)+150.0*(metric[node[i].n_n[0]].xi_z)-25.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_z)+3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_z));
			deter[i].jp = (1.0/256.0)*(3.0*(det[node[node[i].n_n[2]].n_n[2]])-25.0*(det[node[i].n_n[2]])+150.0*(det[i])+150.0*(det[node[i].n_n[0]])-25.0*(det[node[node[i].n_n[0]].n_n[0]])+3.0*(det[node[node[node[i].n_n[0]].n_n[0]].n_n[0]]));
		
			metric[i].zeta_xjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_x)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_x)+150.0*(metric[node[i].n_n[2]].zeta_x)+150.0*(metric[i].zeta_x)-25.0*(metric[node[i].n_n[0]].zeta_x)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_x));
			metric[i].zeta_yjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_y)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_y)+150.0*(metric[node[i].n_n[2]].zeta_y)+150.0*(metric[i].zeta_y)-25.0*(metric[node[i].n_n[0]].zeta_y)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_y));
			metric[i].zeta_zjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_z)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_z)+150.0*(metric[node[i].n_n[2]].zeta_z)+150.0*(metric[i].zeta_z)-25.0*(metric[node[i].n_n[0]].zeta_z)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_z));
			metric[i].eta_xjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].eta_x)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_x)+150.0*(metric[node[i].n_n[2]].eta_x)+150.0*(metric[i].eta_x)-25.0*(metric[node[i].n_n[0]].eta_x)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_x));
			metric[i].eta_yjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].eta_y)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_y)+150.0*(metric[node[i].n_n[2]].eta_y)+150.0*(metric[i].eta_y)-25.0*(metric[node[i].n_n[0]].eta_y)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_y));	
			metric[i].eta_zjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].eta_z)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_z)+150.0*(metric[node[i].n_n[2]].eta_z)+150.0*(metric[i].eta_z)-25.0*(metric[node[i].n_n[0]].eta_z)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_z));	
			metric[i].xi_xjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_x)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_x)+150.0*(metric[node[i].n_n[2]].xi_x)+150.0*(metric[i].xi_x)-25.0*(metric[node[i].n_n[0]].xi_x)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_x));
			metric[i].xi_yjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_y)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_y)+150.0*(metric[node[i].n_n[2]].xi_y)+150.0*(metric[i].xi_y)-25.0*(metric[node[i].n_n[0]].xi_y)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_y));	
			metric[i].xi_zjm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_z)-25.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_z)+150.0*(metric[node[i].n_n[2]].xi_z)+150.0*(metric[i].xi_z)-25.0*(metric[node[i].n_n[0]].xi_z)+3.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_z));	
			deter[i].jm = (1.0/256.0)*(3.0*(det[node[node[node[i].n_n[2]].n_n[2]].n_n[2]])-25.0*(det[node[node[i].n_n[2]].n_n[2]])+150.0*(det[node[i].n_n[2]])+150.0*(det[i])-25.0*(det[node[i].n_n[0]])+3.0*(det[node[node[i].n_n[0]].n_n[0]]));	
			
			metric[i].zeta_xkp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_x)-25.0*(metric[node[i].n_n[5]].zeta_x)+150.0*(metric[i].zeta_x)+150.0*(metric[node[i].n_n[4]].zeta_x)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_x)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_x));
			metric[i].zeta_ykp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_y)-25.0*(metric[node[i].n_n[5]].zeta_y)+150.0*(metric[i].zeta_y)+150.0*(metric[node[i].n_n[4]].zeta_y)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_y)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_y));
			metric[i].zeta_zkp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_z)-25.0*(metric[node[i].n_n[5]].zeta_z)+150.0*(metric[i].zeta_z)+150.0*(metric[node[i].n_n[4]].zeta_z)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_z)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_z));
			metric[i].eta_xkp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_x)-25.0*(metric[node[i].n_n[5]].eta_x)+150.0*(metric[i].eta_x)+150.0*(metric[node[i].n_n[4]].eta_x)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_x)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_x));
			metric[i].eta_ykp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_y)-25.0*(metric[node[i].n_n[5]].eta_y)+150.0*(metric[i].eta_y)+150.0*(metric[node[i].n_n[4]].eta_y)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_y)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_y));
			metric[i].eta_zkp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_z)-25.0*(metric[node[i].n_n[5]].eta_z)+150.0*(metric[i].eta_z)+150.0*(metric[node[i].n_n[4]].eta_z)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_z)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_z));
			metric[i].xi_xkp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_x)-25.0*(metric[node[i].n_n[5]].xi_x)+150.0*(metric[i].xi_x)+150.0*(metric[node[i].n_n[4]].xi_x)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_x)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].xi_x));
			metric[i].xi_ykp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_y)-25.0*(metric[node[i].n_n[5]].xi_y)+150.0*(metric[i].xi_y)+150.0*(metric[node[i].n_n[4]].xi_y)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_y)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].xi_y));
			metric[i].xi_zkp = (1.0/256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_z)-25.0*(metric[node[i].n_n[5]].xi_z)+150.0*(metric[i].xi_z)+150.0*(metric[node[i].n_n[4]].xi_z)-25.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_z)+3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].xi_z));
			deter[i].kp = (1.0/256.0)*(3.0*(det[node[node[i].n_n[5]].n_n[5]])-25.0*(det[node[i].n_n[5]])+150.0*(det[i])+150.0*(det[node[i].n_n[4]])-25.0*(det[node[node[i].n_n[4]].n_n[4]])+3.0*(det[node[node[node[i].n_n[4]].n_n[4]].n_n[4]]));
		
			metric[i].zeta_xkm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_x)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_x)+150.0*(metric[node[i].n_n[5]].zeta_x)+150.0*(metric[i].zeta_x)-25.0*(metric[node[i].n_n[4]].zeta_x)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_x));
			metric[i].zeta_ykm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_y)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_y)+150.0*(metric[node[i].n_n[5]].zeta_y)+150.0*(metric[i].zeta_y)-25.0*(metric[node[i].n_n[4]].zeta_y)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_y));
			metric[i].zeta_zkm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_z)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_z)+150.0*(metric[node[i].n_n[5]].zeta_z)+150.0*(metric[i].zeta_z)-25.0*(metric[node[i].n_n[4]].zeta_z)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_z));
			metric[i].eta_xkm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_x)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_x)+150.0*(metric[node[i].n_n[5]].eta_x)+150.0*(metric[i].eta_x)-25.0*(metric[node[i].n_n[4]].eta_x)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_x));
			metric[i].eta_ykm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_y)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_y)+150.0*(metric[node[i].n_n[5]].eta_y)+150.0*(metric[i].eta_y)-25.0*(metric[node[i].n_n[4]].eta_y)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_y));	
			metric[i].eta_zkm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_z)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_z)+150.0*(metric[node[i].n_n[5]].eta_z)+150.0*(metric[i].eta_z)-25.0*(metric[node[i].n_n[4]].eta_z)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_z));	
			metric[i].xi_xkm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].xi_x)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_x)+150.0*(metric[node[i].n_n[5]].xi_x)+150.0*(metric[i].xi_x)-25.0*(metric[node[i].n_n[4]].xi_x)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_x));
			metric[i].xi_ykm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].xi_y)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_y)+150.0*(metric[node[i].n_n[5]].xi_y)+150.0*(metric[i].xi_y)-25.0*(metric[node[i].n_n[4]].xi_y)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_y));	
			metric[i].xi_zkm = (1.0/256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].xi_z)-25.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_z)+150.0*(metric[node[i].n_n[5]].xi_z)+150.0*(metric[i].xi_z)-25.0*(metric[node[i].n_n[4]].xi_z)+3.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_z));	
			deter[i].km = (1.0/256.0)*(3.0*(det[node[node[node[i].n_n[5]].n_n[5]].n_n[5]])-25.0*(det[node[node[i].n_n[5]].n_n[5]])+150.0*(det[node[i].n_n[5]])+150.0*(det[i])-25.0*(det[node[i].n_n[4]])+3.0*(det[node[node[i].n_n[4]].n_n[4]]));	

		}
		
		/********* Repeating for evaluating metrics at half points for ghost points***********************************/
	//	g_node=sd_node+sd_gh+1;
		for(i=0; i<all_bou_node; i++)
		{			
			gar = 1;
			if (node[all_boundary_nodes[i]].loc > 14 && node[all_boundary_nodes[i]].loc < 27)
			{
				gar++;
			}
			
			if (node[all_boundary_nodes[i]].corner_ID == 0)
			{
				temp1 = 0;
				temp2 = 0;
				temp3 = 0;
				temp4 = 0;
				temp5 = 0;
				temp6 = 0;
				for (temp7 = 0; temp7 < gar; temp7++)
				{
					if (node[node[all_boundary_nodes[i]].n_n[3]].loc == 100 && temp1 == 0)
					{
						j = 3; /********no node on WEST**************/
						k = 1;
						m = 4;	
						gar5 = 1;	
						gar4 = 4;					
						temp1++;
					}
				
					if (node[node[all_boundary_nodes[i]].n_n[2]].loc == 100 && temp2 == 0)
					{
						j = 2; /********no element on SOUTH**************/
						k = 0;
						m = 3;
						gar5 = 2;			
						gar4 = 3;
						temp2++;
					}
					
					if (node[node[all_boundary_nodes[i]].n_n[1]].loc == 100 && temp3 == 0)
					{
						j = 1; /********no element on EAST**************/
						k = 3;
						m = 2;
						gar5 = 1;
						gar4 = 2;
						temp3++;
					}
					
					if (node[node[all_boundary_nodes[i]].n_n[0]].loc == 100 && temp4 == 0)
					{
						j = 0; /********no element on NORTH**************/
						k = 2;
						m = 1;		
						gar5 = 2;	
						gar4 = 1;
						temp4++;
					}
					
					if (node[node[all_boundary_nodes[i]].n_n[4]].loc == 100 && temp5 == 0)
					{
						j = 4; /********no element on NORTH**************/
						k = 5;
						m = 1;		
						gar5 = 3;	
						gar4 = 5;	
						temp5++;
					}
					
					if (node[node[all_boundary_nodes[i]].n_n[5]].loc == 100 && temp6 ==0)
					{
						j = 5; /********no element on NORTH**************/
						k = 4;
						m = 1;		
						gar5 = 3;	
						gar4 = 6;
						temp6++;
					}				
					
					if (node[all_boundary_nodes[i]].loc >= 27 && node[all_boundary_nodes[i]].loc <= 32)
					{
						gar5 = 4;
						gar4 = 2;
					}
					
					switch (gar5)
					{
						case 1:
							//node[all_boundary_nodes[i]].n_n[j] = g_node;
							//node[g_node].n_n[k] = all_boundary_nodes[i];										
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
							deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].im;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
							deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].ip;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
							deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jp;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
							deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jm;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
							deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
							deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;
							
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;						
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;
							
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;						
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
						break;
			
						case 2:
							//node[all_boundary_nodes[i]].n_n[j] = g_node;
							//node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
							deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
							deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;				
							deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jm;				
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
							deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jp;		

							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
							deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
							deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;
					
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;			
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;			
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
						
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;
						
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;			
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;			
								
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
						break;
						
						case 3:
							//node[all_boundary_nodes[i]].n_n[j] = g_node;
							//node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
							deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
							deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;				
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;				
							deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jp;				
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
							deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jm;		

							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
							deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].km;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
							deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;			
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;			
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;			
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;			
								
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;			
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;			
								
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
						break;
						
						case 4:
							for (temp = 0; temp<3; temp++)
							{
								if (temp == 0 && node[i].loc == 28)
								{
									j = 0; /********no element on NORTH**************/
									k = 2;
									m = 1;	
								}
								if (temp == 1 && node[i].loc == 28)
								{
									j = 2; /********no element on SOUTH**************/
									k = 0;
									m = 3;
								}
								if (temp == 0 && node[i].loc == 27)
								{
									j = 4; /********no element on NORTH**************/
									k = 5;
									m = 1;	
								}
								if (temp == 1 && node[i].loc == 27)
								{
									j = 5; /********no element on SOUTH**************/
									k = 4;
									m = 3;
								}
								if (temp == 2)
								{
									j = 1; /********no element on EAST**************/
									k = 3;
									m = 2;
								}
								if (temp <= 1)
								{
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
									deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;
								
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
									deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].im;
							
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
									deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;
								
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
									deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;
									
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
									deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;
								
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
									metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
									metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
									metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
									deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].km;

									
									
									
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
									deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
							
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
									deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
							
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
									deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;
								
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
									deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
									
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
									deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
							
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
									metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
									deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
							
							
									
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
									deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
								
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
									deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
								
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
									deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;
							
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
									deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
									
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
									deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
								
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
									metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
									deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

								}
								if (temp == 2)
								{	
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
									deter[singular[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].im;
							
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
									deter[singular[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;
							
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
									deter[singular[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;
								
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
									deter[singular[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;

									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
									deter[singular[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;
								
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
									metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
									deter[singular[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].km;

								
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
									deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
							
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
									deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
							
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
									deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
								
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
									deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;
									
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
									deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
								
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
									metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
									deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;
			
									
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
									deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
								
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
									deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
									deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
								
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
									deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;
									
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
									deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
								
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
									metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
									deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
								}														
							}						
						break;
					}
					
					switch (gar4)
					{
						case 1:
							metric[all_boundary_nodes[i]].zeta_xjp = metric[all_boundary_nodes[i]].zeta_xjm;
							metric[all_boundary_nodes[i]].eta_xjp = metric[all_boundary_nodes[i]].eta_xjm;
							metric[all_boundary_nodes[i]].xi_xjp = metric[all_boundary_nodes[i]].xi_xjm;
							metric[all_boundary_nodes[i]].zeta_yjp = metric[all_boundary_nodes[i]].zeta_yjm;
							metric[all_boundary_nodes[i]].eta_yjp = metric[all_boundary_nodes[i]].eta_yjm;
							metric[all_boundary_nodes[i]].xi_yjp = metric[all_boundary_nodes[i]].xi_yjm;
							metric[all_boundary_nodes[i]].zeta_zjp = metric[all_boundary_nodes[i]].zeta_zjm;
							metric[all_boundary_nodes[i]].eta_zjp = metric[all_boundary_nodes[i]].eta_zjm;
							metric[all_boundary_nodes[i]].xi_zjp = metric[all_boundary_nodes[i]].xi_zjm;
							deter[all_boundary_nodes[i]].jp = deter[all_boundary_nodes[i]].jm;
						break;
						case 2:
							metric[all_boundary_nodes[i]].zeta_xip = metric[all_boundary_nodes[i]].zeta_xim;
							metric[all_boundary_nodes[i]].eta_xip = metric[all_boundary_nodes[i]].eta_xim;
							metric[all_boundary_nodes[i]].xi_xip = metric[all_boundary_nodes[i]].xi_xim;
							metric[all_boundary_nodes[i]].zeta_yip = metric[all_boundary_nodes[i]].zeta_yim;
							metric[all_boundary_nodes[i]].eta_yip = metric[all_boundary_nodes[i]].eta_yim;
							metric[all_boundary_nodes[i]].xi_yip = metric[all_boundary_nodes[i]].xi_yim;
							metric[all_boundary_nodes[i]].zeta_zip = metric[all_boundary_nodes[i]].zeta_zim;
							metric[all_boundary_nodes[i]].eta_zip = metric[all_boundary_nodes[i]].eta_zim;
							metric[all_boundary_nodes[i]].xi_zip = metric[all_boundary_nodes[i]].xi_zim;
							deter[all_boundary_nodes[i]].ip = deter[all_boundary_nodes[i]].im;
						break;
						case 3:
							metric[all_boundary_nodes[i]].zeta_xjm = metric[all_boundary_nodes[i]].zeta_xjp;
							metric[all_boundary_nodes[i]].eta_xjm = metric[all_boundary_nodes[i]].eta_xjp;
							metric[all_boundary_nodes[i]].xi_xjm = metric[all_boundary_nodes[i]].xi_xjp;
							metric[all_boundary_nodes[i]].zeta_yjm = metric[all_boundary_nodes[i]].zeta_yjp;
							metric[all_boundary_nodes[i]].eta_yjm = metric[all_boundary_nodes[i]].eta_yjp;
							metric[all_boundary_nodes[i]].xi_yjm = metric[all_boundary_nodes[i]].xi_yjp;
							metric[all_boundary_nodes[i]].zeta_zjm = metric[all_boundary_nodes[i]].zeta_zjp;
							metric[all_boundary_nodes[i]].eta_zjm = metric[all_boundary_nodes[i]].eta_zjp;
							metric[all_boundary_nodes[i]].xi_zjm = metric[all_boundary_nodes[i]].xi_zjp;
							deter[all_boundary_nodes[i]].jm = deter[all_boundary_nodes[i]].jp;
						break;
						case 4:
							metric[all_boundary_nodes[i]].zeta_xim = metric[all_boundary_nodes[i]].zeta_xip;
							metric[all_boundary_nodes[i]].eta_xim = metric[all_boundary_nodes[i]].eta_xip;
							metric[all_boundary_nodes[i]].xi_xim = metric[all_boundary_nodes[i]].xi_xip;
							metric[all_boundary_nodes[i]].zeta_yim = metric[all_boundary_nodes[i]].zeta_yip;
							metric[all_boundary_nodes[i]].eta_yim = metric[all_boundary_nodes[i]].eta_yip;
							metric[all_boundary_nodes[i]].xi_yim = metric[all_boundary_nodes[i]].xi_yip;
							metric[all_boundary_nodes[i]].zeta_zim = metric[all_boundary_nodes[i]].zeta_zip;
							metric[all_boundary_nodes[i]].eta_zim = metric[all_boundary_nodes[i]].eta_zip;
							metric[all_boundary_nodes[i]].xi_zim = metric[all_boundary_nodes[i]].xi_zip;
							deter[all_boundary_nodes[i]].im = deter[all_boundary_nodes[i]].ip;
						break;
						case 5:
							metric[all_boundary_nodes[i]].zeta_xkp = metric[all_boundary_nodes[i]].zeta_xkm;
							metric[all_boundary_nodes[i]].eta_xkp = metric[all_boundary_nodes[i]].eta_xkm;
							metric[all_boundary_nodes[i]].xi_xkp = metric[all_boundary_nodes[i]].xi_xkm;
							metric[all_boundary_nodes[i]].zeta_ykp = metric[all_boundary_nodes[i]].zeta_ykm;
							metric[all_boundary_nodes[i]].eta_ykp = metric[all_boundary_nodes[i]].eta_ykm;
							metric[all_boundary_nodes[i]].xi_ykp = metric[all_boundary_nodes[i]].xi_ykm;
							metric[all_boundary_nodes[i]].zeta_zkp = metric[all_boundary_nodes[i]].zeta_zkm;
							metric[all_boundary_nodes[i]].eta_zkp = metric[all_boundary_nodes[i]].eta_zkm;
							metric[all_boundary_nodes[i]].xi_zkp = metric[all_boundary_nodes[i]].xi_zkm;
							deter[all_boundary_nodes[i]].kp = deter[all_boundary_nodes[i]].km;
						break;					
						case 6:
							metric[all_boundary_nodes[i]].zeta_xkm = metric[all_boundary_nodes[i]].zeta_xkp;
							metric[all_boundary_nodes[i]].eta_xkm = metric[all_boundary_nodes[i]].eta_xkp;
							metric[all_boundary_nodes[i]].xi_xkm = metric[all_boundary_nodes[i]].xi_xkp;
							metric[all_boundary_nodes[i]].zeta_ykm = metric[all_boundary_nodes[i]].zeta_ykp;
							metric[all_boundary_nodes[i]].eta_ykm = metric[all_boundary_nodes[i]].eta_ykp;
							metric[all_boundary_nodes[i]].xi_ykm = metric[all_boundary_nodes[i]].xi_ykp;
							metric[all_boundary_nodes[i]].zeta_zkm = metric[all_boundary_nodes[i]].zeta_zkp;
							metric[all_boundary_nodes[i]].eta_zkm = metric[all_boundary_nodes[i]].eta_zkp;
							metric[all_boundary_nodes[i]].xi_zkm = metric[all_boundary_nodes[i]].xi_zkp;
							deter[all_boundary_nodes[i]].km = deter[all_boundary_nodes[i]].kp;
						break;	
					}		
				
				
				}	
			}
			
			else if (node[all_boundary_nodes[i]].corner_ID != 0)
			{
				for (gar=0; gar<3; gar++)
				{
					if (node[all_boundary_nodes[i]].corner_ID == 1)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 2)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 3)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 4)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;	
						}
					}
					
					if (node[all_boundary_nodes[i]].corner_ID == 5)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 6)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 7)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;	
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 8)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;	
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;	
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;	
						}
					}
					switch (gar5)
					{
						case 1:
							//node[all_boundary_nodes[i]].n_n[j] = g_node;
							//node[g_node].n_n[k] = all_boundary_nodes[i];										
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
							deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].im;
												
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
							deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].ip;
													
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
							deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jp;
												
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
							deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jm;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
							deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;
												
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
							deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;

								
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;						
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;

							
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;						
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

						break;
			
						case 2:
							//node[all_boundary_nodes[i]].n_n[j] = g_node;
							//node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
							deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
							deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;				
							deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jm;				
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
							deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jp;		

							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
							deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
							deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;

								
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;			
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;			
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;
		
					
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;			
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;			
								
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
	
						break;
						
						case 3:
							//node[all_boundary_nodes[i]].n_n[j] = g_node;
							//node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
							deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
							deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;				
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;				
							deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jm;				
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
							deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jp;		

							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
							deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;
							
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
							metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
							metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
							deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;
				
						
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;			
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;			
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;			
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;
								
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
							metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
							deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;
		
					
							//node[g_node-1].n_n[j] = g_node;
							//node[g_node].n_n[k] = g_node-1;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;			
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;			
								
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;
							
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;			
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;			
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;			
								
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
							metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
							deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
			
						break;
					}	
				}
			}		
		}
		
		for(i=0; i<all_bou_node; i++)
		{			
			if (node[all_boundary_nodes[i]].loc >= 27 && node[all_boundary_nodes[i]].loc <= 32)
			{
				metric[all_boundary_nodes[i]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[2]].zeta_xjp;
				metric[all_boundary_nodes[i]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[2]].eta_xjp;
				metric[all_boundary_nodes[i]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[2]].xi_xjp;
				metric[all_boundary_nodes[i]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[2]].zeta_yjp;
				metric[all_boundary_nodes[i]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[2]].eta_yjp;
				metric[all_boundary_nodes[i]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[2]].xi_yjp;
				metric[all_boundary_nodes[i]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[2]].zeta_zjp;
				metric[all_boundary_nodes[i]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[2]].eta_zjp;
				metric[all_boundary_nodes[i]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[2]].xi_zjp;
				deter[all_boundary_nodes[i]].jp = 1.0/((metric[all_boundary_nodes[i]].zeta_xjp)*(metric[all_boundary_nodes[i]].eta_yjp)-(metric[all_boundary_nodes[i]].eta_xjp)*(metric[all_boundary_nodes[i]].zeta_yjp));
								
				metric[all_boundary_nodes[i]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[3]].zeta_xip;
				metric[all_boundary_nodes[i]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[3]].eta_xip;
				metric[all_boundary_nodes[i]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[3]].xi_xip;
				metric[all_boundary_nodes[i]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[3]].zeta_yip;
				metric[all_boundary_nodes[i]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[3]].eta_yip;
				metric[all_boundary_nodes[i]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[3]].xi_yip;
				metric[all_boundary_nodes[i]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[3]].zeta_zip;
				metric[all_boundary_nodes[i]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[3]].eta_zip;
				metric[all_boundary_nodes[i]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[3]].xi_zip;
				deter[all_boundary_nodes[i]].ip = 1.0/((metric[all_boundary_nodes[i]].zeta_xip)*(metric[all_boundary_nodes[i]].eta_yip)-(metric[all_boundary_nodes[i]].eta_xip)*(metric[all_boundary_nodes[i]].zeta_yip));
						
				metric[all_boundary_nodes[i]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[0]].zeta_xjm;
				metric[all_boundary_nodes[i]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[0]].eta_xjm;
				metric[all_boundary_nodes[i]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[0]].xi_xjm;
				metric[all_boundary_nodes[i]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[0]].zeta_yjm;
				metric[all_boundary_nodes[i]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[0]].eta_yjm;
				metric[all_boundary_nodes[i]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[0]].xi_yjm;
				metric[all_boundary_nodes[i]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[0]].zeta_zjm;
				metric[all_boundary_nodes[i]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[0]].eta_zjm;
				metric[all_boundary_nodes[i]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[0]].xi_zjm;
				deter[all_boundary_nodes[i]].jm = 1.0/((metric[all_boundary_nodes[i]].zeta_xjm)*(metric[all_boundary_nodes[i]].eta_yjm)-(metric[all_boundary_nodes[i]].eta_xjm)*(metric[all_boundary_nodes[i]].zeta_yjm));
			
			}
		}
		
		
		/*******************************************************malloc for variables********************************************/

		u = (double **)malloc(5*sizeof(double *));
		v = (double **)malloc(5*sizeof(double *));
		w = (double **)malloc(5*sizeof(double *));
		rho = (double **)malloc(5*sizeof(double *));
		p = (double **)malloc(5*sizeof(double *));
		t = (double **)malloc(5*sizeof(double *));
		mu = (double **)malloc(5*sizeof(double *));
		a = (double **)malloc(5*sizeof(double *));
		e = (double **)malloc(5*sizeof(double *));
	
		tauzz = (double *)malloc(nodes*sizeof(double));
		tauee = (double *)malloc(nodes*sizeof(double));
		tauxx = (double *)malloc(nodes*sizeof(double));
		
		tauze = (double *)malloc(nodes*sizeof(double));
		tauzx = (double *)malloc(nodes*sizeof(double));
		tauex = (double *)malloc(nodes*sizeof(double));
		
		TAU_SGS_XX = (double *)malloc(nodes*sizeof(double));
		TAU_SGS_YY = (double *)malloc(nodes*sizeof(double));
		TAU_SGS_ZZ = (double *)malloc(nodes*sizeof(double));
		TAU_SGS_XY = (double *)malloc(nodes*sizeof(double));
		TAU_SGS_XZ = (double *)malloc(nodes*sizeof(double));
		TAU_SGS_YZ = (double *)malloc(nodes*sizeof(double));
		H_SGS_X = (double *)malloc(nodes*sizeof(double));
		H_SGS_Y = (double *)malloc(nodes*sizeof(double));
		H_SGS_Z = (double *)malloc(nodes*sizeof(double));
		D_SGS_X = (double *)malloc(nodes*sizeof(double));
		D_SGS_Y = (double *)malloc(nodes*sizeof(double));
		D_SGS_Z = (double *)malloc(nodes*sizeof(double));
		DUCROS = (double *)malloc(nodes*sizeof(double));
		
		qz = (double *)malloc(nodes*sizeof(double));
		qe = (double *)malloc(nodes*sizeof(double));
		qx = (double *)malloc(nodes*sizeof(double));
		
		del_cfl = (double *)malloc(nodes*sizeof(double));
		v_dash1 = (double *)malloc(nodes*sizeof(double));
		
		for (i=0; i<5; i++)
		{
			u[i] = (double *)malloc(nodes*sizeof(double));
			v[i] = (double *)malloc(nodes*sizeof(double));
			w[i] = (double *)malloc(nodes*sizeof(double));
			rho[i] = (double *)malloc(nodes*sizeof(double));
			p[i] = (double *)malloc(nodes*sizeof(double));
			t[i] = (double *)malloc(nodes*sizeof(double));
			mu[i] = (double *)malloc(nodes*sizeof(double));
			a[i] = (double *)malloc(nodes*sizeof(double));
			e[i] = (double *)malloc(nodes*sizeof(double));
		}

		/*******************************malloc for solver variables********************************************/	
		roe_rho_ip = (double *)malloc(nodes*sizeof(double));
		roe_u_ip = (double *)malloc(nodes*sizeof(double));
		roe_v_ip = (double *)malloc(nodes*sizeof(double));
		roe_w_ip = (double *)malloc(nodes*sizeof(double));
		roe_h_ip = (double *)malloc(nodes*sizeof(double));
		
		final_U = (double *)malloc(5*sizeof(double));
		
		roe_rho_im = (double *)malloc(nodes*sizeof(double));
		roe_u_im = (double *)malloc(nodes*sizeof(double));
		roe_v_im = (double *)malloc(nodes*sizeof(double));
		roe_w_im = (double *)malloc(nodes*sizeof(double));
		roe_h_im = (double *)malloc(nodes*sizeof(double));
			
		roe_rho_jp = (double *)malloc(nodes*sizeof(double));
		roe_u_jp = (double *)malloc(nodes*sizeof(double));
		roe_v_jp = (double *)malloc(nodes*sizeof(double));
		roe_w_jp = (double *)malloc(nodes*sizeof(double));
		roe_h_jp = (double *)malloc(nodes*sizeof(double));
		
		roe_rho_jm = (double *)malloc(nodes*sizeof(double));
		roe_u_jm = (double *)malloc(nodes*sizeof(double));
		roe_v_jm = (double *)malloc(nodes*sizeof(double));
		roe_w_jm = (double *)malloc(nodes*sizeof(double));
		roe_h_jm = (double *)malloc(nodes*sizeof(double));
		
		roe_rho_kp = (double *)malloc(nodes*sizeof(double));
		roe_u_kp = (double *)malloc(nodes*sizeof(double));
		roe_v_kp = (double *)malloc(nodes*sizeof(double));
		roe_w_kp = (double *)malloc(nodes*sizeof(double));
		roe_h_kp = (double *)malloc(nodes*sizeof(double));
		
		roe_rho_km = (double *)malloc(nodes*sizeof(double));
		roe_u_km = (double *)malloc(nodes*sizeof(double));
		roe_v_km = (double *)malloc(nodes*sizeof(double));
		roe_w_km = (double *)malloc(nodes*sizeof(double));
		roe_h_km = (double *)malloc(nodes*sizeof(double));
		
		roe_R = (double *)malloc(nodes*sizeof(double));
		roe_a = (double *)malloc(nodes*sizeof(double));

		div =  (RESIDUAL **)malloc(nodes*sizeof(RESIDUAL *));

		Qi_iminus = (double **)malloc(nodes*sizeof(double *));
		Qi_iplus = (double **)malloc(nodes*sizeof(double *));
		Qj_iminus = (double **)malloc(nodes*sizeof(double *));
		Qj_iplus = (double **)malloc(nodes*sizeof(double *));
		Qk_iminus = (double **)malloc(nodes*sizeof(double *));
		Qk_iplus = (double **)malloc(nodes*sizeof(double *));
		
		Qi_iminus_n = (double **)malloc(nodes*sizeof(double *));
		Qi_iplus_n = (double **)malloc(nodes*sizeof(double *));
		Qj_iminus_n = (double **)malloc(nodes*sizeof(double *));
		Qj_iplus_n = (double **)malloc(nodes*sizeof(double *));
		Qk_iminus_n = (double **)malloc(nodes*sizeof(double *));
		Qk_iplus_n = (double **)malloc(nodes*sizeof(double *));

		w_Qip = (double **)malloc(5*sizeof(double *));
		w_Qinp = (double **)malloc(5*sizeof(double *));
		w_Qim = (double **)malloc(5*sizeof(double *));
		w_Qinm = (double **)malloc(5*sizeof(double *));
		w_Qjp = (double **)malloc(5*sizeof(double *));
		w_Qjnp = (double **)malloc(5*sizeof(double *));
		w_Qjm = (double **)malloc(5*sizeof(double *));
		w_Qjnm = (double **)malloc(5*sizeof(double *));
		w_Qkp = (double **)malloc(5*sizeof(double *));
		w_Qknp = (double **)malloc(5*sizeof(double *));
		w_Qkm = (double **)malloc(5*sizeof(double *));
		w_Qknm = (double **)malloc(5*sizeof(double *));
		
		
		W_Qip = (double **)malloc(5*sizeof(double *));
		W_Qinp = (double **)malloc(5*sizeof(double *));
		W_Qim = (double **)malloc(5*sizeof(double *));
		W_Qinm = (double **)malloc(5*sizeof(double *));
		W_Qjp = (double **)malloc(5*sizeof(double *));
		W_Qjnp = (double **)malloc(5*sizeof(double *));
		W_Qjm = (double **)malloc(5*sizeof(double *));
		W_Qjnm = (double **)malloc(5*sizeof(double *));
		W_Qkp = (double **)malloc(5*sizeof(double *));
		W_Qknp = (double **)malloc(5*sizeof(double *));
		W_Qkm = (double **)malloc(5*sizeof(double *));
		W_Qknm = (double **)malloc(5*sizeof(double *));
		
		L = (double ***)malloc(nodes*sizeof(double **));
		U = (double ***)malloc(nodes*sizeof(double **));
		E = (double **)malloc(nodes*sizeof(double *));
		F = (double **)malloc(nodes*sizeof(double *));
		G = (double **)malloc(nodes*sizeof(double *));
		Ev = (double **)malloc(nodes*sizeof(double *));
		Fv = (double **)malloc(nodes*sizeof(double *));
		Gv = (double **)malloc(nodes*sizeof(double *));
		E1 = (double **)malloc(nodes*sizeof(double *));
		F1 = (double **)malloc(nodes*sizeof(double *));
		G1 = (double **)malloc(nodes*sizeof(double *));
		Ev1 = (double **)malloc(nodes*sizeof(double *));
		Fv1 = (double **)malloc(nodes*sizeof(double *));
		Gv1 = (double **)malloc(nodes*sizeof(double *));
		dF = (double *)malloc(5*sizeof(double *));
		dE = (double *)malloc(5*sizeof(double *));
		dG = (double *)malloc(5*sizeof(double *));
		dFv = (double *)malloc(5*sizeof(double *));
		dEv = (double *)malloc(5*sizeof(double *));
		dGv = (double *)malloc(5*sizeof(double *));
		
		
		Qi_half_m = (double **)malloc(5*sizeof(double * ));
		Qi_halfn_m = (double **)malloc(5*sizeof(double * ));
		Qi_half_p = (double **)malloc(5*sizeof(double * ));
		Qi_half_np = (double **)malloc(5*sizeof(double * ));
		Qj_half_m = (double **)malloc(5*sizeof(double * ));
		Qj_halfn_m = (double **)malloc(5*sizeof(double * ));
		Qj_half_p = (double **)malloc(5*sizeof(double * ));
		Qj_half_np = (double **)malloc(5*sizeof(double * ));
		Qk_half_m = (double **)malloc(5*sizeof(double * ));
		Qk_halfn_m = (double **)malloc(5*sizeof(double * ));
		Qk_half_p = (double **)malloc(5*sizeof(double * ));
		Qk_half_np = (double **)malloc(5*sizeof(double * ));
		
		IS_Qim = (double **)malloc(5*sizeof(double *));
		IS_Qinm = (double **)malloc(5*sizeof(double *));
		IS_Qip = (double **)malloc(5*sizeof(double *));
		IS_Qinp = (double **)malloc(5*sizeof(double *));
		IS_Qjm = (double **)malloc(5*sizeof(double *));
		IS_Qjnm = (double **)malloc(5*sizeof(double *));
		IS_Qjp = (double **)malloc(5*sizeof(double *));
		IS_Qjnp = (double **)malloc(5*sizeof(double *));
		IS_Qkm = (double **)malloc(5*sizeof(double *));
		IS_Qknm = (double **)malloc(5*sizeof(double *));
		IS_Qkp = (double **)malloc(5*sizeof(double *));
		IS_Qknp = (double **)malloc(5*sizeof(double *));
		
		r_eigen_Qip = (double**)malloc(5*sizeof(double*));
		r_eigen_Qjp = (double**)malloc(5*sizeof(double*));
		r_eigen_Qkp = (double**)malloc(5*sizeof(double*));
		l_eigen_Qip = (double**)malloc(5*sizeof(double*));
		l_eigen_Qjp = (double**)malloc(5*sizeof(double*));
		l_eigen_Qkp = (double**)malloc(5*sizeof(double*));
		
		r_eigen_Qim = (double**)malloc(5*sizeof(double*));
		r_eigen_Qjm = (double**)malloc(5*sizeof(double*));
		r_eigen_Qkm = (double**)malloc(5*sizeof(double*));
		l_eigen_Qim = (double**)malloc(5*sizeof(double*));
		l_eigen_Qjm = (double**)malloc(5*sizeof(double*));
		l_eigen_Qkm = (double**)malloc(5*sizeof(double*));
		
		for (i=0; i<nodes; i++)
		{
			
			Qi_iminus[i] = (double *)malloc(5*sizeof(double));
			Qi_iplus[i] = (double *)malloc(5*sizeof(double));
			Qj_iminus[i] = (double *)malloc(5*sizeof(double));
			Qj_iplus[i] = (double *)malloc(5*sizeof(double));
			Qk_iminus[i] = (double *)malloc(5*sizeof(double));
			Qk_iplus[i] = (double *)malloc(5*sizeof(double));

			div[i] =  (RESIDUAL*)malloc(10*sizeof(RESIDUAL));
			
			Qi_iminus_n[i] = (double *)malloc(5*sizeof(double));
			Qi_iplus_n[i] = (double *)malloc(5*sizeof(double));
			Qj_iminus_n[i] = (double *)malloc(5*sizeof(double));
			Qj_iplus_n[i] = (double *)malloc(5*sizeof(double));
			Qk_iminus_n[i] = (double *)malloc(5*sizeof(double));
			Qk_iplus_n[i] = (double *)malloc(5*sizeof(double));
				
			L[i] = (double**)malloc(5*sizeof(double*));
			U[i] = (double**)malloc(5*sizeof(double*));
			E[i] = (double*)malloc(5*sizeof(double));
			F[i] = (double*)malloc(5*sizeof(double));
			G[i] = (double*)malloc(5*sizeof(double));
			Ev[i] = (double*)malloc(5*sizeof(double));
			Fv[i] = (double*)malloc(5*sizeof(double));
			Gv[i] = (double*)malloc(5*sizeof(double));
			E1[i] = (double*)malloc(5*sizeof(double));
			F1[i] = (double*)malloc(5*sizeof(double));
			G1[i] = (double*)malloc(5*sizeof(double));
			Ev1[i] = (double*)malloc(5*sizeof(double));
			Fv1[i] = (double*)malloc(5*sizeof(double));
			Gv1[i] = (double*)malloc(5*sizeof(double));
		}
		
		for (i=0;i<5; i++)
		{
			w_Qip[i] = (double*)malloc(5*sizeof(double));
			w_Qinp[i] = (double*)malloc(5*sizeof(double));
			w_Qim[i] = (double*)malloc(5*sizeof(double));
			w_Qinm[i] = (double*)malloc(5*sizeof(double));
			w_Qjp[i] = (double*)malloc(5*sizeof(double));
			w_Qjnp[i] = (double*)malloc(5*sizeof(double));
			w_Qjm[i] = (double*)malloc(5*sizeof(double));
			w_Qjnm[i] = (double*)malloc(5*sizeof(double));
			w_Qkp[i] = (double*)malloc(5*sizeof(double));
			w_Qknp[i] = (double*)malloc(5*sizeof(double));
			w_Qkm[i] = (double*)malloc(5*sizeof(double));
			w_Qknm[i] = (double*)malloc(5*sizeof(double));
		
			W_Qip[i] = (double*)malloc(5*sizeof(double));
			W_Qinp[i] = (double*)malloc(5*sizeof(double));
			W_Qim[i] = (double*)malloc(5*sizeof(double));
			W_Qinm[i] = (double*)malloc(5*sizeof(double));
			W_Qjp[i] = (double*)malloc(5*sizeof(double));
			W_Qjnp[i] = (double*)malloc(5*sizeof(double));
			W_Qjm[i] = (double*)malloc(5*sizeof(double));
			W_Qjnm[i] = (double*)malloc(5*sizeof(double));
			W_Qkp[i] = (double*)malloc(5*sizeof(double));
			W_Qknp[i] = (double*)malloc(5*sizeof(double));
			W_Qkm[i] = (double*)malloc(5*sizeof(double));
			W_Qknm[i] = (double*)malloc(5*sizeof(double));
			
			IS_Qim[i] = (double*)malloc(5*sizeof(double));
			IS_Qinm[i] = (double*)malloc(5*sizeof(double));
			IS_Qip[i] = (double*)malloc(5*sizeof(double));
			IS_Qinp[i] = (double*)malloc(5*sizeof(double));
			IS_Qjm[i] = (double*)malloc(5*sizeof(double));
			IS_Qjnm[i] = (double*)malloc(5*sizeof(double));
			IS_Qjp[i] = (double*)malloc(5*sizeof(double));
			IS_Qjnp[i] = (double*)malloc(5*sizeof(double));
			IS_Qkm[i] = (double*)malloc(5*sizeof(double));
			IS_Qknm[i] = (double*)malloc(5*sizeof(double));
			IS_Qkp[i] = (double*)malloc(5*sizeof(double));
			IS_Qknp[i] = (double*)malloc(5*sizeof(double));
			Qi_half_m[i] = (double*)malloc(5*sizeof(double));
			Qi_halfn_m[i] = (double*)malloc(5*sizeof(double));
			Qi_half_p[i] = (double*)malloc(5*sizeof(double));
			Qi_half_np[i] = (double*)malloc(5*sizeof(double));
			Qj_half_m[i] = (double*)malloc(5*sizeof(double));
			Qj_halfn_m[i] = (double*)malloc(5*sizeof(double));
			Qj_half_p[i] = (double*)malloc(5*sizeof(double));
			Qj_half_np[i] = (double*)malloc(5*sizeof(double));
			Qk_half_m[i] = (double*)malloc(5*sizeof(double));
			Qk_halfn_m[i] = (double*)malloc(5*sizeof(double));
			Qk_half_p[i] = (double*)malloc(5*sizeof(double));
			Qk_half_np[i] = (double*)malloc(5*sizeof(double));
			
			r_eigen_Qip[i] = (double*)malloc(5*sizeof(double));
			r_eigen_Qjp[i] = (double*)malloc(5*sizeof(double));
			r_eigen_Qkp[i] = (double*)malloc(5*sizeof(double));
			l_eigen_Qip[i] = (double*)malloc(5*sizeof(double));
			l_eigen_Qjp[i] = (double*)malloc(5*sizeof(double));
			l_eigen_Qkp[i] = (double*)malloc(5*sizeof(double));
			l_eigen_Qim[i] = (double*)malloc(5*sizeof(double));
			l_eigen_Qjm[i] = (double*)malloc(5*sizeof(double));
			l_eigen_Qkm[i] = (double*)malloc(5*sizeof(double));
			
			r_eigen_Qim[i] = (double*)malloc(5*sizeof(double));
			r_eigen_Qjm[i] = (double*)malloc(5*sizeof(double));
			r_eigen_Qkm[i] = (double*)malloc(5*sizeof(double));
		}

		for (i=0; i<nodes;i++)
		{
			for(j=0; j<5; j++)
			{
				L[i][j] = (double*)malloc(5*sizeof(double));
				U[i][j] = (double*)malloc(5*sizeof(double));
			}	
		}
		
		int *temp_a;
		temp_a = (int*)malloc((sd_node+10)*sizeof(int));
		/******************************Reading restart file**********************************************/
		sprintf(filename,"restart_file_%d.neu",restart_num);
		fnode = fopen(filename,"rt");
		if (fnode != NULL )
		{
			i=1;	
			while(fgets(line, 1000, fnode) != NULL)
			{
				sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %d", &u[0][d_node[i].local], &v[0][d_node[i].local], &w[0][d_node[i].local], &rho[0][d_node[i].local], &p[0][d_node[i].local], &t[0][d_node[i].local], &e[0][d_node[i].local], &mu[0][d_node[i].local], &itera);
				i++;
				if(i > NUMNP)
				{
					break;
				}
			}	
			itera++;
			itera++;
			fclose(fnode);		
		}
		else
		{
			itera = 1;
		}
		
		free(d_node);
		d_node = NULL;
		/***************************SENDING LOCAL NODE DATA TO MASTER************************************/
		MPI_Send(&sd_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);
		
		position = 0;
		MPI_Pack_size(sd_node, MPI_INT, MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		for (i=1; i<=sd_node; i++)
		{
			temp_a[i] = node[i].global;
			MPI_Pack(&temp_a[i],1 ,MPI_INT ,buffer ,memsize ,&position ,MPI_COMM_WORLD);
		}
		MPI_Send(buffer, memsize, MPI_PACKED,0,myrank,MPI_COMM_WORLD);
		free (buffer);
		buffer = NULL;
		//free(temp_a);
		
		/*************************************************************************************************
									Flow variables initializing
		*************************************************************************************************/
	//	Reyl = 7500.0;
	//	int iterations = 1000000;
		int back_p;
		j = 0;
		initial = 0;
		if (itera != 1)
		{
			initial = 1;
			back_p = (itera-1180000)/100000;
			back_pressure = back_pressure+0.1*back_p;
		}
		intialise(u, v, w, rho, p, t, mu, e, g_node, inlet_node, outlet_node, wall_node, boundary_node, CD, wal_node, initial, jacobian, metric, j, a, \
		nodes, NUMNP, out_node, inl_node, node, bou_node, all_bou_node, all_boundary_nodes, div, sd_node, sd_wal_node, sd_bou_node, sd_out_node, sd_inl_node, sd_wall_node, \
		sd_boundary_node, sd_outlet_node, sd_inlet_node, singular);
	
		temp1 = sd_node;
		position = 0;
		for(i=0; i<neigh_pro; i++)
		{
			MPI_Pack_size(recv_c[c[i]]*9,MPI_DOUBLE,MPI_COMM_WORLD,&memsize1);
			buffer2 = malloc(memsize1);                                          /***********carefull with buffer1 and buffer2******************/
			MPI_Pack_size(proc_node[c[i]]*9,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);
			position = 0;
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&u[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&v[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}	
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&w[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&p[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&t[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&rho[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&e[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for(j=0; j<recv_c[c[i]];j++)
			{
				MPI_Pack(&mu[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
				
			MPI_Sendrecv(buffer2, memsize1, MPI_PACKED, c[i], c[i], buffer, memsize, MPI_PACKED, c[i], myrank, MPI_COMM_WORLD, &status);
			position = 0;
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&u[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&v[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&w[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&p[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&t[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&rho[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&e[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			for(j=0; j<proc_node[c[i]];j++)
			{
				MPI_Unpack(buffer,memsize,&position,&mu[0][recv_b[i][j]],1,MPI_DOUBLE,MPI_COMM_WORLD);
			}
			free(buffer);
			buffer = NULL;
			free(buffer2);
			buffer2 = NULL;
		}
		
		j =0;
		initial = 1;
		intialise(u, v, w, rho, p, t, mu, e, g_node, inlet_node, outlet_node, wall_node, boundary_node, CD, wal_node, initial, jacobian, metric, j, a, \
		nodes, NUMNP, out_node, inl_node, node, bou_node, all_bou_node, all_boundary_nodes, div, sd_node, sd_wal_node, sd_bou_node, sd_out_node, sd_inl_node, sd_wall_node, \
		sd_boundary_node, sd_outlet_node, sd_inlet_node, singular);	
		
		sd_node= temp1;
		for (k=1; k<=no_of_tip_send; k++)
		{
			u[0][tip[k].node] = 0.0;
			v[0][tip[k].node] = 0.0;
			w[0][tip[k].node] = 0.0;
			p[0][tip[k].node] = p[0][node[tip[k].node].n_n[3]];
			t[0][tip[k].node] = t[0][node[tip[k].node].n_n[3]];
			rho[0][tip[k].node] = (1.4*Mach*Mach)*(p[0][tip[k].node]/t[0][tip[k].node]);
			e[0][tip[k].node] = p[0][tip[k].node]/(0.4*rho[0][tip[k].node]);
		}
			
		for (k=1; k<=no_of_tip_recv; k++)
		{
			u[0][tip_recv[k].node] = 0.0;
			v[0][tip_recv[k].node] = 0.0;
			w[0][tip_recv[k].node] = 0.0;
			p[0][tip_recv[k].node] = p[0][node[tip_recv[k].node].n_n[3]];
			t[0][tip_recv[k].node] = t[0][node[tip_recv[k].node].n_n[3]];
			rho[0][tip_recv[k].node] = (1.4*Mach*Mach)*(p[0][tip_recv[k].node]/t[0][tip_recv[k].node]);
			e[0][tip_recv[k].node] = p[0][tip_recv[k].node]/(0.4*rho[0][tip_recv[k].node]);
		}
		
		j = 0;			
		sd_node= temp1;

		if (itera == 1)
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
			
			MPI_Send(buffer, memsize, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);			
			free(buffer);
			buffer = NULL;	
		}
		
		MPI_Barrier(comm_slaves);
		
		/***********************************************************************************************************************************
																	SOLVER LOOP
		************************************************************************************************************************************/	
		solver(CD, jacobian, metric, U, F1, E1, G1, Fv1, Ev1, Gv1, E, F, G, Ev, Fv, Gv, roe_rho_ip, roe_u_ip, roe_v_ip, roe_w_ip, roe_h_ip, roe_rho_im, roe_u_im, roe_v_im, \
		roe_w_im, roe_h_im, roe_rho_jp, roe_u_jp, roe_v_jp, roe_w_jp, roe_h_jp, roe_rho_jm, roe_u_jm, roe_v_jm, roe_w_jm, roe_h_jm, roe_rho_kp, roe_u_kp, roe_v_kp, roe_w_kp, \
		roe_h_kp, roe_rho_km, roe_u_km, roe_v_km, roe_w_km, roe_h_km, r_eigen_Qip, r_eigen_Qjp, r_eigen_Qkp, r_eigen_Qim, r_eigen_Qjm, r_eigen_Qkm, l_eigen_Qip, l_eigen_Qjp, \
		l_eigen_Qkp, Qi_iplus, Qi_iminus, Qj_iplus, Qj_iminus, Qk_iplus, Qk_iminus, Qi_half_p, Qi_half_np, Qi_half_m, Qi_halfn_m, Qj_half_p, Qj_half_np, Qj_half_m, \
		Qj_halfn_m, Qk_half_p, Qk_half_np, Qk_half_m, Qk_halfn_m, IS_Qim, IS_Qinm, IS_Qjm, IS_Qjnm,  IS_Qkm, IS_Qknm, IS_Qip, IS_Qinp, IS_Qjp, IS_Qjnp, IS_Qkp, IS_Qknp, w_Qip, \
		w_Qinp, w_Qim, w_Qinm, w_Qjp, w_Qjnp, w_Qjm, w_Qjnm, w_Qkp, w_Qknp, w_Qkm, w_Qknm, W_Qip, W_Qinp, W_Qim, W_Qinm, W_Qjp, W_Qjnp, W_Qjm, W_Qjnm, W_Qkp, W_Qknp, W_Qkm, W_Qknm, \
		Qi_iplus_half, Qi_iminus_half, Qj_iplus_half, Qj_iminus_half,  Qk_iplus_half, Qk_iminus_half, dF, dE, dG, dFv, dEv, dGv, L, NUMNP, NELEM, rho, u, v, w, e, p, t, \
		a, mu, g_node, g_elem, cor, all_bou_node, inlet, outlet, wall, boundary, iterations, tauzz, tauee,  tauxx, tauze, tauzx, tauex, qz, qe, qx, det, all_boundary_nodes, Reyl, \
		Qi_iplus_n, Qi_iminus_n, Qj_iplus_n, Qj_iminus_n, Qk_iplus_n, Qk_iminus_n, final_U, inl, out, wal, bou, inlet_node, outlet_node, wall_node, boundary_node, inl_elem, \
		out_elem, bou_elem, wal_elem, nodes, wal_node, node, roe_R, roe_a, out_node, inl_node, bou_node, v_dash1, del_cfl, div, l_eigen_Qim, l_eigen_Qjm, l_eigen_Qkm, sd_node, loc_dat, \
		recv_b, c, proc_node, neigh_pro, recv_c, itera, no_of_tip_send, no_of_tip_recv, tip, tip_recv, sd_wal_node, sd_bou_node, sd_out_node, sd_inl_node, sd_wall_node, \
		sd_boundary_node, sd_outlet_node, sd_inlet_node, singular, Q, deter, TAU_SGS_XY, TAU_SGS_YZ, TAU_SGS_XZ, H_SGS_X, H_SGS_Y, H_SGS_Z, D_SGS_X, D_SGS_Y, D_SGS_Z, TAU_SGS_XX, TAU_SGS_YY, TAU_SGS_ZZ, DUCROS);
	}
	//MPI_Barrier(MPI_COMM_WORLD);

	/**********************DATA GATHERING AND FILE WRITING AND PRINTING AFTER ITERATIONS************************************************/
	if (myrank == 0)
	{
	/*	MPI_Pack_size(NUMNP*4, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Pack(&node[i].x, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].y, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].z, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		}		
		for (i=1; i<=size-1; i++)
		{
			MPI_Send(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD);	
		}
		free(buffer);
		buffer = NULL;	
		
		
		MPI_Pack_size(NUMNP*19, MPI_INT, MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Pack(&node[i].n_n[0], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[1], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[2], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[3], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[4], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[5], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[0], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[1], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[2], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[3], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[4], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[5], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[6], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].e[7], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].loc, 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].val, 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].ID, 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].corner_ID, 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
		}		
		for (i=1; i<=size-1; i++)
		{
			MPI_Send(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD);	
		}
		free(buffer);
		buffer = NULL;	
		*/
		//system("rm -rf proc* boundary* nodes_proc* outlet* wall* inlet* ");
		double temp_u, temp_v, temp_w, temp_p, temp_t, temp_e, temp_mu, temp_rho;
		double *max_div_u, *max_div_v, *max_div_e;
		int *sd_count, **b, glob, iter, *sd_gh, temp_nodes;
		
		//printf("\n%d\n",_LINE_);
		
		
		printf("\nfreeing memory of temp arrays\n");
	/*	
		free(jacobian);
		jacobian = NULL;
		free(metric);
		metric = NULL;
		free(det);
		det = NULL;
		free(deter);
		deter = NULL;
		*/
		max_div_u = (double *)malloc(size*sizeof(double));
		max_div_v = (double *)malloc(size*sizeof(double));
		max_div_e = (double *)malloc(size*sizeof(double));
		sd_count = (int*)malloc(size*sizeof(int));
		b = (int **)malloc(size*sizeof(int *));
		sd_gh = (int *)malloc(size*sizeof(int));
		
		u = (double **)malloc(5*sizeof(double *));
		v = (double **)malloc(5*sizeof(double *));
		w = (double **)malloc(5*sizeof(double *));
		rho = (double **)malloc(5*sizeof(double *));
		p = (double **)malloc(5*sizeof(double *));
		t = (double **)malloc(5*sizeof(double *));
		mu = (double **)malloc(5*sizeof(double *));
		e = (double **)malloc(5*sizeof(double *));		
		DUCROS = (double *)malloc((NUMNP+10)*sizeof(double));
		
		for (i=0; i<5; i++)
		{
			u[i] = (double *)malloc(nodes*sizeof(double));
			v[i] = (double *)malloc(nodes*sizeof(double));
			w[i] = (double *)malloc(nodes*sizeof(double));
			rho[i] = (double *)malloc(nodes*sizeof(double));
			p[i] = (double *)malloc(nodes*sizeof(double));
			t[i] = (double *)malloc(nodes*sizeof(double));
			mu[i] = (double *)malloc(nodes*sizeof(double));
			//a[i] = (double *)malloc(nodes*sizeof(double));
			e[i] = (double *)malloc(nodes*sizeof(double));
			//ki[i] = (double *)malloc(nodes*sizeof(double));
			
		}
		
		/*****************************************************************************************************************/
		int *neigh_pro, **all_c_list;
		neigh_pro = (int *)malloc((size+1)*sizeof(int));
		all_c_list = (int **)malloc((size+1)*sizeof(int *));
		
		for (i=1; i<=size-1; i++)
		{			
			MPI_Recv(&neigh_pro[i],1, MPI_INT, i, i, MPI_COMM_WORLD,&status);		
			
			all_c_list[i] = (int *)malloc((neigh_pro[i]+1)*sizeof(int));
			position = 0;
			MPI_Pack_size(neigh_pro[i], MPI_INT,MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);
			MPI_Recv(buffer,memsize,MPI_PACKED,i,i,MPI_COMM_WORLD,&status);		
			for (j=1; j<=neigh_pro[i]; j++)
			{
				MPI_Unpack(buffer,memsize,&position,&glob,1,MPI_INT,MPI_COMM_WORLD);
				all_c_list[i][j] = glob;
			}
			free(buffer);
			buffer = NULL;
		}

		
		
		for (i=1; i<=size-1; i++)
		{	
			for (k=1; k<=size-1; k++)
			{
				MPI_Send(&neigh_pro[k], 1, MPI_INT, i, i, MPI_COMM_WORLD);
				position = 0;
				MPI_Pack_size(neigh_pro[k], MPI_INT,MPI_COMM_WORLD,&memsize);
				buffer = malloc(memsize);
				
				for (j=1; j<=neigh_pro[k]; j++)
				{
					MPI_Pack(&all_c_list[k][j],1,MPI_INT,buffer,memsize,&position,MPI_COMM_WORLD);
				}
				
				MPI_Send(buffer,memsize,MPI_PACKED,i,i,MPI_COMM_WORLD);				
				free(buffer);
				buffer = NULL;
			}
		}
		
		/********************************************************************************************************************/	
	
		for (i=1; i<=size-1; i++)
		{
			MPI_Recv(&sd_gh[i],1, MPI_INT, i, i, MPI_COMM_WORLD,&status);		
			position = 0;
			MPI_Pack_size(sd_gh[i], MPI_INT,MPI_COMM_WORLD,&memsize);
			buffer = malloc(memsize);
			MPI_Recv(buffer,memsize,MPI_PACKED,i,i,MPI_COMM_WORLD,&status);		
			for (j=1; j< sd_gh[i]; j++)
			{
				MPI_Unpack(buffer,memsize,&position,&glob,1,MPI_INT,MPI_COMM_WORLD);
				node[glob].req = i;
			}
			free(buffer);
			buffer = NULL;
		}
		
	/*	
		MPI_Pack_size(NUMNP*9, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
		buffer = malloc(memsize);
		position = 0;
		for(i=1; i<=NUMNP; i++)
		{
			MPI_Pack(&jacobian[i].x_zeta, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].x_eta, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].x_xi, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_zeta, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_eta, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_xi, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_zeta, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_eta, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_xi, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			
		}
		for (i=1; i<=size-1; i++)
		{
			MPI_Send(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD);	
		}
		free(buffer);
		buffer = NULL;	
		*/
		
		 free(jacobian);
		jacobian = NULL;
		free(metric);
		metric = NULL;
		free(det);
		det = NULL;
		free(deter);
		deter = NULL;
		
		
		//printf("\n7378\n");
		
		FILE *fp;
		sprintf(filename,"proc_request.plt");
		fp= fopen(filename,"w");
		fprintf(fp,"TITLE = \"Node file\"\n");
		fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"Req\", \"Proc\", \"loc\",\n");
		fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",NUMNP,NELEM);
		for(i=0;i<NUMNP;i++)
		{
			fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\n", node[i+1].x, node[i+1].y, node[i+1].z, node[i+1].req, node[i+1].proc, node[i+1].loc);
		}
		fprintf(fp,"\n\n\n");
	
		for(i=1; i<=NELEM; i++)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", CD[i].connect[0], CD[i].connect[1], CD[i].connect[5], CD[i].connect[4], CD[i].connect[3], CD[i].connect[2], CD[i].connect[6], CD[i].connect[7] );
		}			
		fclose(fp);
		
		for (i=1; i<=size-1; i++)
		{	
			MPI_Recv(&sd_count[i],1, MPI_INT, i, i, MPI_COMM_WORLD,&status);
			position = 0;
			MPI_Pack_size(sd_count[i], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = malloc(memsize);
			MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);	
			b[i] = (int *)malloc((sd_count[i]+10)*sizeof(int));		
			for (j=1; j<=sd_count[i]; j++)
			{
				MPI_Unpack(buffer,memsize,&position,&glob,1,MPI_INT,MPI_COMM_WORLD);
				b[i][j] = glob;
			}
			free(buffer);
			buffer = NULL;
		}
		
		sprintf(filename,"restart_file_%d.neu",restart_num);
		fp = fopen(filename,"rt");
		if (fp != NULL)
		{
			fgets(line, 1000, fp);
			sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %d", &temp_u, &temp_v, &temp_w, &temp_p, &temp_t, &temp_rho, &temp_e, &temp_mu, &itera);
			itera++;
			fclose (fp);
		}	
		else
		{
			itera = 0;
		}
		
		printf("deleting all temperory files\n");
		system("rm -rf *.txt *.o node_neighbour.neu elem_neighbour.neu");
		printf("calculations started\n");
		
		restart_num = 0;
	//	iterations = 2000000;
		for (iter=itera; iter<iterations; iter++)
		{	
		//	restart_num++;
			//printf("\n7434\n");
			if (iter > itera)
			{
			//	printf("\n7437\n");
				for (i=1; i<=size-1; i++)
				{	
					position =0;
					MPI_Pack_size(3, MPI_DOUBLE, MPI_COMM_WORLD,&memsize);
					buffer = malloc(memsize);
					MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);	
					MPI_Unpack(buffer, memsize, &position, &temp_d, 1, MPI_DOUBLE, MPI_COMM_WORLD);		
					max_div_u[i] = temp_d;
					MPI_Unpack(buffer, memsize, &position, &temp_d, 1, MPI_DOUBLE, MPI_COMM_WORLD);		
					max_div_v[i] = temp_d;	
					MPI_Unpack(buffer, memsize, &position, &temp_d, 1, MPI_DOUBLE, MPI_COMM_WORLD);		
					max_div_e[i] = temp_d;
				
					free(buffer);
					buffer = NULL;
				}
				
			//	printf("\n7453\n");
				
				for (i=1; i<= size-1; i++)
				{
					for (i=1; i<=size-1; i++)
					{
						if (max_div_u[i] < max_div_u[i+1])
						{
							temp_d = max_div_u[i+1];
							max_div_u[i+1] = max_div_u[i];
							max_div_u[i] = temp_d;
						}
						if (max_div_v[i] < max_div_v[i+1])
						{
							temp_d = max_div_v[i+1];
							max_div_v[i+1] = max_div_v[i];
							max_div_v[i] = temp_d;
						}
						if (max_div_e[i] < max_div_e[i+1])
						{
							temp_d = max_div_e[i+1];
							max_div_e[i+1] = max_div_e[i];
							max_div_e[i] = temp_d;
						}
					}
				}
				printf("iteration = %d   max_div_u = %e max_div_v = %e max_div_e = %e \n", iter-1, max_div_u[1], max_div_v[1], max_div_e[1]);
			}

			if (iter%100 == 0 )
			{				
				restart_num++;
				for (i=1; i<=size-1; i++)
				{
					position = 0;
					MPI_Pack_size(sd_count[i]*10.0,MPI_DOUBLE,MPI_COMM_WORLD,&memsize);
					buffer = malloc(memsize);
					temp = sd_count[i]*8.0;
				
					MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);	
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &u[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);	
					}
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &v[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &w[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &rho[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);	
					}
					for (j=1; j<=sd_count[i]; j++)
					{		
						MPI_Unpack(buffer, memsize, &position, &p[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}	
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &t[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);	
					}
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &e[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);	
					}
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &mu[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j=1; j<=sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &DUCROS[b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					free(buffer);
					buffer = NULL;
				}	
				if (iter%1000 == 0 )
				{
					writepltfile(NELEM, NUMNP, node, CD, u, v, w, rho, p, t, e, mu, iter, DUCROS);
					
					/* if (restart_num > 2)
					{
						restart_num = 1;
					}	 */			
					//restart_file(NELEM, NUMNP, node, CD, u, v, w, rho, p, t, e, mu, iter, restart_num);
					//writefile(NELEM, NUMNP, node, CD, tauzz, tauze, tauee, qz, qe, u_zeta, v_zeta, t_zeta, u_eta, v_eta, t_eta);
					
				 	sprintf(filename,"nohup ./preplot nodefile_%d.dat > preplot.out",iter);
					system(filename);
					system("rm -rf preplot.out *.dat"); 
				}
				
				if (restart_num > 2)
				{
					restart_num = 1;
				}
				
				restart_file(NELEM, NUMNP, node, CD, u, v, w, rho, p, t, e, mu, iter, restart_num);
			}			
		}		
	}	
}
