#include <stdio.h>
#include<stdlib.h> 
#include<math.h>
#include"function.h"

void intialise(double **u, double **v, double **w, double **rho, double **p, double **t, double **mu, double **e, int g_node, int *inlet_node, int *outlet_node, int *wall_node, \
int *boundary_node, ELEM *CD, int wal_node, int initial, JACOB *jacobian, TRANSFORMED *metric,int j, double **a, int nodes, int NUMNP, int out_node, int inl_node, MNODE *node, \
int bou_node, int all_bou_node, int *all_boundary_nodes, RESIDUAL **div, int sd_node, int sd_wal_node,int sd_bou_node, int sd_out_node, int sd_inl_node, int *sd_wall_node, \
int *sd_boundary_node, int *sd_outlet_node, int *sd_inlet_node, singul *singular)
{
	int i,l,k,h;
	double Cd, Cdr1, Cdr2, Cdr3, Cdc1, Cdc2, Cdc3, determinant, angle, angle2, Vt, Vn, Vt2, Vn2, Vn_sign;
	Cd = 0.07;
	Cdr1 = 0.0123; 
	Cdr2 = 0.0279;
	Cdr3 = 0.0264;
	Cdc1 = 0.0121;
	Cdc2 = 0.0138;
	Cdc3 = 0.0;
	
	if (initial == 0)
	{
		for (i=1; i<=sd_node; i++)
		{				
			u[j][i] = 1.0;
			v[j][i] = 0.0;
			w[j][i] = 0.0;
			p[j][i] = 1.0/(1.4*Mach*Mach);
			t[j][i] = 1.0;
			rho[j][i] = (1.4*Mach*Mach)*(p[j][i]/t[j][i]);
			e[j][i] = p[j][i]/(0.4*rho[j][i]);
			a[j][i] = sqrt(1.4*p[j][i]/rho[j][i]);																
		}		
	}
	
/****************************************************OUTLET********************************************************/		
	for (i=0; i<out_node; i++)
	{	
		if (node[outlet_node[i]].loc > 0 && node[outlet_node[i]].loc <= 6)
		{	
			if (node[outlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[outlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[outlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[outlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[outlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[outlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[j][outlet_node[i]]= u[j][node[outlet_node[i]].n_n[k]];	
			v[j][outlet_node[i]]= v[j][node[outlet_node[i]].n_n[k]];	
			w[j][outlet_node[i]]= w[j][node[outlet_node[i]].n_n[k]];	
			if (node[outlet_node[i]].y < 0.0 || node[outlet_node[i]].y > 380.0)
			{
				p[j][outlet_node[i]]= p[j][node[outlet_node[i]].n_n[k]];	
			}
			if (node[outlet_node[i]].y > 0.0 && node[outlet_node[i]].y < 380.0)
			{
				p[j][outlet_node[i]]= back_pressure*(1.0/(1.4*Mach*Mach));	
			}
			
			t[j][outlet_node[i]]= t[j][node[outlet_node[i]].n_n[k]];	
			rho[j][outlet_node[i]]= (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
			a[j][outlet_node[i]] = sqrt(1.4*p[j][outlet_node[i]]/rho[j][outlet_node[i]]); 

			u[j][node[outlet_node[i]].n_n[l]] = u[j][outlet_node[i]];
			v[j][node[outlet_node[i]].n_n[l]] = v[j][outlet_node[i]];
			w[j][node[outlet_node[i]].n_n[l]] = w[j][outlet_node[i]];
			p[j][node[outlet_node[i]].n_n[l]] = p[j][outlet_node[i]];
			rho[j][node[outlet_node[i]].n_n[l]] = rho[j][outlet_node[i]];
			t[j][node[outlet_node[i]].n_n[l]] = t[j][outlet_node[i]];	
			a[j][node[outlet_node[i]].n_n[l]] = a[j][outlet_node[i]];
			e[j][node[outlet_node[i]].n_n[l]] = e[j][outlet_node[i]];
			
			u[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = u[j][outlet_node[i]];
			v[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = v[j][outlet_node[i]];
			w[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = w[j][outlet_node[i]];
			p[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = p[j][outlet_node[i]];
			rho[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = rho[j][outlet_node[i]];
			t[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = t[j][outlet_node[i]];	
			a[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = a[j][outlet_node[i]];
			e[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = e[j][outlet_node[i]];
	
			u[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][outlet_node[i]];
			v[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][outlet_node[i]];
			w[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][outlet_node[i]];
			p[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][outlet_node[i]];
			rho[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][outlet_node[i]];
			t[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][outlet_node[i]];	
			a[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][outlet_node[i]];
			e[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][outlet_node[i]];
			
			 
		/* 	u[j][node[outlet_node[i]].n_n[l]] = 2.0*u[j][outlet_node[i]]-u[j][node[outlet_node[i]].n_n[k]];
			v[j][node[outlet_node[i]].n_n[l]] = 2.0*v[j][outlet_node[i]]-v[j][node[outlet_node[i]].n_n[k]];
			w[j][node[outlet_node[i]].n_n[l]] = 2.0*w[j][outlet_node[i]]-w[j][node[outlet_node[i]].n_n[k]];
			p[j][node[outlet_node[i]].n_n[l]] = 2.0*p[j][outlet_node[i]]-p[j][node[outlet_node[i]].n_n[k]];
			rho[j][node[outlet_node[i]].n_n[l]] = 2.0*rho[j][outlet_node[i]]-rho[j][node[outlet_node[i]].n_n[k]];
			t[j][node[outlet_node[i]].n_n[l]] = 2.0*t[j][outlet_node[i]]-t[j][node[outlet_node[i]].n_n[k]];	
			a[j][node[outlet_node[i]].n_n[l]] = a[j][outlet_node[i]];
			e[j][node[outlet_node[i]].n_n[l]] = e[j][outlet_node[i]];
			
			u[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*u[j][node[outlet_node[i]].n_n[l]]-u[j][outlet_node[i]];
			v[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*v[j][node[outlet_node[i]].n_n[l]]-v[j][outlet_node[i]];
			w[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*w[j][node[outlet_node[i]].n_n[l]]-w[j][outlet_node[i]];
			p[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*p[j][node[outlet_node[i]].n_n[l]]-p[j][outlet_node[i]];
			rho[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*rho[j][node[outlet_node[i]].n_n[l]]-rho[j][outlet_node[i]];
			t[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*t[j][node[outlet_node[i]].n_n[l]]-t[j][outlet_node[i]];
			a[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = a[j][outlet_node[i]];
			e[j][node[node[outlet_node[i]].n_n[l]].n_n[l]] = e[j][outlet_node[i]];
	
			u[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*u[j][node[node[outlet_node[i]].n_n[l]].n_n[l]]-u[j][node[outlet_node[i]].n_n[l]];
			v[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*v[j][node[node[outlet_node[i]].n_n[l]].n_n[l]]-v[j][node[outlet_node[i]].n_n[l]];
			w[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*w[j][node[node[outlet_node[i]].n_n[l]].n_n[l]]-w[j][node[outlet_node[i]].n_n[l]];
			p[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*p[j][node[node[outlet_node[i]].n_n[l]].n_n[l]]-p[j][node[outlet_node[i]].n_n[l]];
			rho[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*rho[j][node[node[outlet_node[i]].n_n[l]].n_n[l]]-rho[j][node[outlet_node[i]].n_n[l]];
			t[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*t[j][node[node[outlet_node[i]].n_n[l]].n_n[l]]-t[j][node[outlet_node[i]].n_n[l]];
			a[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][outlet_node[i]];
			e[j][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][outlet_node[i]];
			 */
		}
		else if (node[outlet_node[i]].loc == 7 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);

		}
		else if (node[outlet_node[i]].loc == 8 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 9 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 10 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		
		else if (node[outlet_node[i]].loc == 11 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 12 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 13 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 14 )
		{
			u[j][outlet_node[i]] = u[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);

		}
		else if (node[outlet_node[i]].loc == 15 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 16 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 17 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 18 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		
		else if (node[outlet_node[i]].loc == 19 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 20 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 21 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 22 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 23 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 24 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 25 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 26 )
		{
			u[j][outlet_node[i]] = u[j][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			v[j][outlet_node[i]] = v[j][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			w[j][outlet_node[i]] = w[j][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			p[j][outlet_node[i]] = p[j][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			t[j][outlet_node[i]] = t[j][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			rho[j][outlet_node[i]] = (1.4*Mach*Mach)*(p[j][outlet_node[i]]/t[j][outlet_node[i]]);
			e[j][outlet_node[i]] = p[j][outlet_node[i]]/(0.4*rho[j][outlet_node[i]]);
		}
	}

/****************************************************INLET********************************************************/
	for (i=0; i< inl_node; i++)
	{	
		if (node[inlet_node[i]].loc > 0 && node[inlet_node[i]].loc <= 6)
		{	
			if (node[inlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[inlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[inlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[inlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[inlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[inlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}
		
			u[j][inlet_node[i]]= 1.0;		
			v[j][inlet_node[i]]= 0.0;	
			w[j][inlet_node[i]]= 0.0;				
			p[j][inlet_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[j][inlet_node[i]]= 1.0;		
			rho[j][inlet_node[i]]= (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
			a[j][inlet_node[i]] = sqrt(1.4*p[j][inlet_node[i]]/rho[j][inlet_node[i]]); 

			u[j][node[inlet_node[i]].n_n[l]] = u[j][inlet_node[i]];
			v[j][node[inlet_node[i]].n_n[l]] = v[j][inlet_node[i]];
			w[j][node[inlet_node[i]].n_n[l]] = w[j][inlet_node[i]];
			p[j][node[inlet_node[i]].n_n[l]] = p[j][inlet_node[i]];
			rho[j][node[inlet_node[i]].n_n[l]] = rho[j][inlet_node[i]];
			t[j][node[inlet_node[i]].n_n[l]] = t[j][inlet_node[i]];	
			a[j][node[inlet_node[i]].n_n[l]] = a[j][inlet_node[i]];
			e[j][node[inlet_node[i]].n_n[l]] = e[j][inlet_node[i]];
			//mu[j][node[inlet_node[i]].n_n[l]] = mu[j][inlet_node[i]];
			
			u[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = u[j][inlet_node[i]];
			v[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = v[j][inlet_node[i]];
			w[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = w[j][inlet_node[i]];
			p[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = p[j][inlet_node[i]];
			rho[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = rho[j][inlet_node[i]];
			t[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = t[j][inlet_node[i]];	
			a[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = a[j][inlet_node[i]];
			e[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = e[j][inlet_node[i]];
			//mu[j][node[node[inlet_node[i]].n_n[l]].n_n[l]] = mu[j][inlet_node[i]];
	
			u[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][inlet_node[i]];
			v[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][inlet_node[i]];
			w[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][inlet_node[i]];
			p[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][inlet_node[i]];
			rho[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][inlet_node[i]];
			t[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][inlet_node[i]];	
			a[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][inlet_node[i]];
			e[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][inlet_node[i]];
			//mu[j][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][inlet_node[i]];
		}
		else if (node[inlet_node[i]].loc == 7 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 8 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 9 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 10 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		
		else if (node[inlet_node[i]].loc == 11 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 12 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 13 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 14 )
		{
			u[j][inlet_node[i]] = u[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 15 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 16 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 17 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 18 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		
		else if (node[inlet_node[i]].loc == 19 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 20 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 21 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 22 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 23 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 24 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 25 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 26 )
		{
			u[j][inlet_node[i]] = u[j][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			v[j][inlet_node[i]] = v[j][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			w[j][inlet_node[i]] = w[j][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			p[j][inlet_node[i]] = p[j][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			t[j][inlet_node[i]] = t[j][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			rho[j][inlet_node[i]] = (1.4*Mach*Mach)*(p[j][inlet_node[i]]/t[j][inlet_node[i]]);
			e[j][inlet_node[i]] = p[j][inlet_node[i]]/(0.4*rho[j][inlet_node[i]]);
		}
	}
	
	/****************************************************BOUNDARY********************************************************/
	for (i=0; i< bou_node; i++)
	{	
		if (node[boundary_node[i]].loc > 0 && node[boundary_node[i]].loc <= 6)
		{	
			if (node[boundary_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[boundary_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[boundary_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[boundary_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[boundary_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[boundary_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

/*			u[j][boundary_node[i]]= 1.0;		
			v[j][boundary_node[i]]= 0.0;	
			w[j][boundary_node[i]]= 0.0;				
			p[j][boundary_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[j][boundary_node[i]]= 1.0;		
			rho[j][boundary_node[i]]= (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
			a[j][boundary_node[i]] = sqrt(1.4*p[j][boundary_node[i]]/rho[j][boundary_node[i]]); 
*/
			u[j][boundary_node[i]] = u[j][node[boundary_node[i]].n_n[k]];
			v[j][boundary_node[i]] = v[j][node[boundary_node[i]].n_n[k]];
			w[j][boundary_node[i]] = w[j][node[boundary_node[i]].n_n[k]];
			p[j][boundary_node[i]] = p[j][node[boundary_node[i]].n_n[k]];
			rho[j][boundary_node[i]] = rho[j][node[boundary_node[i]].n_n[k]];
			t[j][boundary_node[i]] = t[j][node[boundary_node[i]].n_n[k]];	
			a[j][boundary_node[i]] = a[j][node[boundary_node[i]].n_n[k]];
			e[j][boundary_node[i]] = e[j][node[boundary_node[i]].n_n[k]];
			//mu[j][boundary_node[i]] = mu[j][boundary_node[i]];
			
			u[j][node[boundary_node[i]].n_n[l]] = u[j][node[boundary_node[i]].n_n[k]];
			v[j][node[boundary_node[i]].n_n[l]] = v[j][node[boundary_node[i]].n_n[k]];
			w[j][node[boundary_node[i]].n_n[l]] = w[j][node[boundary_node[i]].n_n[k]];
			p[j][node[boundary_node[i]].n_n[l]] = p[j][node[boundary_node[i]].n_n[k]];
			rho[j][node[boundary_node[i]].n_n[l]] = rho[j][node[boundary_node[i]].n_n[k]];
			t[j][node[boundary_node[i]].n_n[l]] = t[j][node[boundary_node[i]].n_n[k]];	
			a[j][node[boundary_node[i]].n_n[l]] = a[j][node[boundary_node[i]].n_n[k]];
			e[j][node[boundary_node[i]].n_n[l]] = e[j][node[boundary_node[i]].n_n[k]];
			//mu[j][node[boundary_node[i]].n_n[l]] = mu[j][boundary_node[i]];
			
			u[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = u[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			v[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = v[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			w[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = w[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			p[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = p[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			rho[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = rho[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			t[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = t[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];	
			a[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = a[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			e[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = e[j][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			//mu[j][node[node[boundary_node[i]].n_n[l]].n_n[l]] = mu[j][boundary_node[i]];
	
			u[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];	
			a[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[j][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][boundary_node[i]];
		}
		else if (node[boundary_node[i]].loc == 7 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);

		}
		else if (node[boundary_node[i]].loc == 8 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 9 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 10 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		
		else if (node[boundary_node[i]].loc == 11 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 12 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 13 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 14 )
		{
			u[j][boundary_node[i]] = u[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);

		}
		else if (node[boundary_node[i]].loc == 15 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 16 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 17 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 18 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		
		else if (node[boundary_node[i]].loc == 19 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 20 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 21 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 22 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 23 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 24 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 25 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 26 )
		{
			u[j][boundary_node[i]] = u[j][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			v[j][boundary_node[i]] = v[j][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			w[j][boundary_node[i]] = w[j][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			p[j][boundary_node[i]] = p[j][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			t[j][boundary_node[i]] = t[j][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			rho[j][boundary_node[i]] = (1.4*Mach*Mach)*(p[j][boundary_node[i]]/t[j][boundary_node[i]]);
			e[j][boundary_node[i]] = p[j][boundary_node[i]]/(0.4*rho[j][boundary_node[i]]);
		}		
	}
	
	/*****************************************************************************************************************************************************************/
	
	for (i=0; i< wal_node; i++)
	{	
		if (node[wall_node[i]].loc > 0 && node[wall_node[i]].loc <= 6)
		{	
			if (node[wall_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[wall_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[wall_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[wall_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[wall_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[wall_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[wall_node[i]].n_n[k]];
			//p[j][wall_node[i]]= 1.0/(1.4*Mach*Mach);
			t[j][wall_node[i]]= t[j][node[wall_node[i]].n_n[k]];	
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);		
			a[j][wall_node[i]] = sqrt(1.4*p[j][wall_node[i]]/rho[j][wall_node[i]]);	
		
			u[j][node[wall_node[i]].n_n[l]] = (-1.0)*u[j][node[wall_node[i]].n_n[k]];
			v[j][node[wall_node[i]].n_n[l]] = (-1.0)*v[j][node[wall_node[i]].n_n[k]];
			w[j][node[wall_node[i]].n_n[l]] = (-1.0)*w[j][node[wall_node[i]].n_n[k]];
			p[j][node[wall_node[i]].n_n[l]] = p[j][node[wall_node[i]].n_n[k]];
			rho[j][node[wall_node[i]].n_n[l]] = rho[j][node[wall_node[i]].n_n[k]];
			t[j][node[wall_node[i]].n_n[l]] = t[j][node[wall_node[i]].n_n[k]];
			a[j][node[wall_node[i]].n_n[l]] = a[j][node[wall_node[i]].n_n[k]];
			e[j][node[wall_node[i]].n_n[l]] = e[j][node[wall_node[i]].n_n[k]];
			//mu[j][node[wall_node[i]].n_n[l]] = mu[j][wall_node[i]];
		
			u[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			v[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			w[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			p[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			rho[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			t[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			a[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			e[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[j][node[node[wall_node[i]].n_n[k]].n_n[k]];
			//mu[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[j][wall_node[i]];
		
			u[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			a[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][wall_node[i]];
			
			/*******************************BLEED*********************************************/
		/* 	if ( (node[wall_node[i]].x >= 1120.0 && node[wall_node[i]].x <= 1220.0 && node[wall_node[i]].y < 354.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1350.0 && node[wall_node[i]].x <= 1410.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1015.238 && node[wall_node[i]].x <= 1250.0 && node[wall_node[i]].loc == 1) || (node[wall_node[i]].x >= 1350.0 && node[wall_node[i]].x <= 1482.0 && node[wall_node[i]].loc == 1))
//			if ((node[wall_node[i]].x >= 1226.5+700.0 && node[wall_node[i]].x <= 1411.224+700.0 && node[wall_node[i]].y < 674.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1411.224+700.0 && node[wall_node[i]].x <= 1564.39+700.0 && node[wall_node[i]].loc == 3))	
			{
			//	 if (node[wall_node[i]].x >= 720.0 && node[wall_node[i]].x <= 1226.5 && node[wall_node[i]].y < 354.0 && node[wall_node[i]].loc == 3)
			//	{
			//		Cd = Cdr1;
			//	} 
				 
				if (node[wall_node[i]].x >= 1120.0 && node[wall_node[i]].x <= 1220.0 && node[wall_node[i]].y < 354.0 && node[wall_node[i]].loc == 3)
				{
					Cd = Cdr2;
				}
				
				if (node[wall_node[i]].x >= 1350.0 && node[wall_node[i]].x <= 1410.0 && node[wall_node[i]].loc == 3)
				{
					Cd = Cdr3;
				}
				
				if (node[wall_node[i]].x >= 1015.238 && node[wall_node[i]].x <= 1250.0 && node[wall_node[i]].loc == 1)
				{
					Cd = Cdc1;
				}
		 	
				if (node[wall_node[i]].x >= 1350.0 && node[wall_node[i]].x <= 1482.0 && node[wall_node[i]].loc == 1)
				{
					Cd = Cdc2;
				}
				
				if (node[wall_node[i]].loc == 3)
				{
					Vn_sign = -1.0;
				}
				
				if (node[wall_node[i]].loc == 1)
				{
					Vn_sign = 1.0;
				}
				
				rho[j][wall_node[i]]= rho[j][node[wall_node[i]].n_n[k]];
				angle = atan((node[node[wall_node[i]].n_n[1]].y-node[wall_node[i]].y)/(node[node[wall_node[i]].n_n[1]].x-node[wall_node[i]].x));
				angle2 = atan((node[node[node[wall_node[i]].n_n[k]].n_n[1]].y-node[node[wall_node[i]].n_n[k]].y)/(node[node[node[wall_node[i]].n_n[k]].n_n[1]].x-node[node[wall_node[i]].n_n[k]].x));
				determinant = cos(angle)*cos(angle)+sin(angle)*sin(angle);

				Vt = u[j][node[wall_node[i]].n_n[k]]*cos(angle)+v[j][node[wall_node[i]].n_n[k]]*sin(angle);
				Vn = Vn_sign*Cd*sqrt(1.4*p[j][node[wall_node[i]].n_n[k]]/rho[j][node[wall_node[i]].n_n[k]]);
				Vn2 = v[j][node[wall_node[i]].n_n[k]]*cos(angle2)-u[j][node[wall_node[i]].n_n[k]]*sin(angle2);
				
				p[j][wall_node[i]] = p[j][node[wall_node[i]].n_n[k]]-2.0*rho[j][wall_node[i]]*Vn*(Vn-Vn2);				
				u[j][wall_node[i]] = (cos(angle)/determinant)*Vt-(sin(angle)/determinant)*Vn;
				v[j][wall_node[i]] = (sin(angle)/determinant)*Vt+(cos(angle)/determinant)*Vn;
				w[j][wall_node[i]]= 0.0;
				t[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/rho[j][wall_node[i]]);
				e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);		
				a[j][wall_node[i]] = sqrt(1.4*p[j][wall_node[i]]/rho[j][wall_node[i]]);	
				
				u[j][node[wall_node[i]].n_n[l]] = u[j][wall_node[i]];
				v[j][node[wall_node[i]].n_n[l]] = v[j][wall_node[i]];
				w[j][node[wall_node[i]].n_n[l]] = w[j][wall_node[i]];
				p[j][node[wall_node[i]].n_n[l]] = p[j][wall_node[i]];
				rho[j][node[wall_node[i]].n_n[l]] = rho[j][wall_node[i]];
				t[j][node[wall_node[i]].n_n[l]] = t[j][wall_node[i]];	
				a[j][node[wall_node[i]].n_n[l]] = a[j][wall_node[i]];
				e[j][node[wall_node[i]].n_n[l]] = e[j][wall_node[i]];
				
				u[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = u[j][wall_node[i]];
				v[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = v[j][wall_node[i]];
				w[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = w[j][wall_node[i]];
				p[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[j][wall_node[i]];
				rho[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[j][wall_node[i]];
				t[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[j][wall_node[i]];	
				a[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[j][wall_node[i]];
				e[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[j][wall_node[i]];
		
				u[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][wall_node[i]];
				v[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][wall_node[i]];
				w[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][wall_node[i]];
				p[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][wall_node[i]];
				rho[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][wall_node[i]];
				t[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][wall_node[i]];	
				a[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][wall_node[i]];
				e[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][wall_node[i]];
				
			} 		
			 */
		}
		else if (node[wall_node[i]].loc == 7 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);

		}
		else if (node[wall_node[i]].loc == 8 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 9 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 10 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		
		else if (node[wall_node[i]].loc == 11 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 12 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 13 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 14 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 15 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[wall_node[i]].n_n[5]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[5]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 16 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[wall_node[i]].n_n[1]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[1]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 17 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]]= p[j][node[node[wall_node[i]].n_n[4]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[4]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 18 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[3]].n_n[0]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[3]].n_n[0]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		
		else if (node[wall_node[i]].loc == 19 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[5]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[5]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 20 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[1]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[1]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 21 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[4]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[4]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 22 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[3]].n_n[2]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[3]].n_n[2]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 23 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[1]].n_n[5]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[1]].n_n[5]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 24 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[1]].n_n[4]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[1]].n_n[4]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 25 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[3]].n_n[4]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[3]].n_n[4]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 26 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[node[wall_node[i]].n_n[3]].n_n[5]];
			t[j][wall_node[i]] = t[j][node[node[wall_node[i]].n_n[3]].n_n[5]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 29 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[singular[wall_node[i]].n_n[3]].n_n[4]];
			t[j][wall_node[i]] = t[j][node[singular[wall_node[i]].n_n[3]].n_n[4]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 30 )
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][node[singular[wall_node[i]].n_n[3]].n_n[5]];
			t[j][wall_node[i]] = t[j][node[singular[wall_node[i]].n_n[3]].n_n[5]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);
		}
		
		if (node[wall_node[i]].ID == 10 && node[wall_node[i]].loc != 29 && node[wall_node[i]].loc != 30)
		{
			u[j][wall_node[i]]= 0.0;
			v[j][wall_node[i]]= 0.0;
			w[j][wall_node[i]]= 0.0;
			p[j][wall_node[i]] = p[j][singular[wall_node[i]].n_n[3]];
			t[j][wall_node[i]] = t[j][singular[wall_node[i]].n_n[3]];
			rho[j][wall_node[i]] = (1.4*Mach*Mach)*(p[j][wall_node[i]]/t[j][wall_node[i]]);
			e[j][wall_node[i]] = p[j][wall_node[i]]/(0.4*rho[j][wall_node[i]]);		
			
			for (h=0;h<3;h++)
			{
				if (h==0)
				{
					l = 2; /********no element on SOUTH**************/
					k = 0;
				}
				if (h==1) 
				{
					l = 0; /********no element on NORTH**************/
					k = 2;
				}
				if (h==2)
				{
					l = 1; /********no element on EAST**************/
					k = 3;										
				}
				
				if (h<=1)
				{
					u[j][node[wall_node[i]].n_n[l]] = (-1.0)*u[j][singular[wall_node[i]].n_n[k]];
					v[j][node[wall_node[i]].n_n[l]] = (-1.0)*v[j][singular[wall_node[i]].n_n[k]];
					w[j][node[wall_node[i]].n_n[l]] = (-1.0)*w[j][singular[wall_node[i]].n_n[k]];
					p[j][node[wall_node[i]].n_n[l]] = p[j][singular[wall_node[i]].n_n[k]];
					rho[j][node[wall_node[i]].n_n[l]] = rho[j][singular[wall_node[i]].n_n[k]];
					t[j][node[wall_node[i]].n_n[l]] = t[j][singular[wall_node[i]].n_n[k]];	
					a[j][node[wall_node[i]].n_n[l]] = a[j][wall_node[i]];
					e[j][node[wall_node[i]].n_n[l]] = e[j][singular[wall_node[i]].n_n[k]];
					//mu[j][node[wall_node[i]].n_n[l]] = mu[j][wall_node[i]];
				
					u[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					v[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					w[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					p[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					rho[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					t[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];	
					a[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[j][wall_node[i]];
					e[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					//mu[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[j][wall_node[i]];
		
					u[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][wall_node[i]];
					e[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][wall_node[i]];	
				}
				if (h == 2)
				{
					u[j][node[wall_node[i]].n_n[l]] = (-1.0)*u[j][singular[wall_node[i]].n_n[k]];
					v[j][node[wall_node[i]].n_n[l]] = (-1.0)*v[j][singular[wall_node[i]].n_n[k]];
					w[j][node[wall_node[i]].n_n[l]] = (-1.0)*w[j][singular[wall_node[i]].n_n[k]];
					p[j][node[wall_node[i]].n_n[l]] = p[j][singular[wall_node[i]].n_n[k]];
					rho[j][node[wall_node[i]].n_n[l]] = rho[j][singular[wall_node[i]].n_n[k]];
					t[j][node[wall_node[i]].n_n[l]] = t[j][singular[wall_node[i]].n_n[k]];	
					a[j][node[wall_node[i]].n_n[l]] = a[j][wall_node[i]];
					e[j][node[wall_node[i]].n_n[l]] = e[j][singular[wall_node[i]].n_n[k]];
					//mu[j][node[wall_node[i]].n_n[l]] = mu[j][wall_node[i]];
				
					u[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					v[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					w[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					p[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					rho[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					t[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];	
					a[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[j][wall_node[i]];
					e[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[j][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					//mu[j][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[j][wall_node[i]];
		
					u[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][wall_node[i]];
					e[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[j][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][wall_node[i]];	
				}
			}
		}
	}
	
	/****************************************************sd_outlet********************************************************/		
	for (i=0; i<sd_out_node; i++)
	{	
		if (node[sd_outlet_node[i]].loc > 0 && node[sd_outlet_node[i]].loc <= 6)
		{	
			if (node[sd_outlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_outlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_outlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_outlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_outlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_outlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[j][sd_outlet_node[i]]= u[j][node[sd_outlet_node[i]].n_n[k]];	
			v[j][sd_outlet_node[i]]= v[j][node[sd_outlet_node[i]].n_n[k]];	
			w[j][sd_outlet_node[i]]= w[j][node[sd_outlet_node[i]].n_n[k]];	
			if (node[sd_outlet_node[i]].y < 0.0 || node[sd_outlet_node[i]].y > 380.0)
			{
				p[j][sd_outlet_node[i]]= p[j][node[sd_outlet_node[i]].n_n[k]];	
			}
			if (node[sd_outlet_node[i]].y > 0.0 && node[sd_outlet_node[i]].y < 380.0)
			{
				p[j][sd_outlet_node[i]]= back_pressure*(1.0/(1.4*Mach*Mach));	
			}
			
			t[j][sd_outlet_node[i]]= t[j][node[sd_outlet_node[i]].n_n[k]];	
			rho[j][sd_outlet_node[i]]= (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
			a[j][sd_outlet_node[i]] = sqrt(1.4*p[j][sd_outlet_node[i]]/rho[j][sd_outlet_node[i]]); 
			
			u[j][node[sd_outlet_node[i]].n_n[l]] = u[j][sd_outlet_node[i]];
			v[j][node[sd_outlet_node[i]].n_n[l]] = v[j][sd_outlet_node[i]];
			w[j][node[sd_outlet_node[i]].n_n[l]] = w[j][sd_outlet_node[i]];
			p[j][node[sd_outlet_node[i]].n_n[l]] = p[j][sd_outlet_node[i]];
			rho[j][node[sd_outlet_node[i]].n_n[l]] = rho[j][sd_outlet_node[i]];
			t[j][node[sd_outlet_node[i]].n_n[l]] = t[j][sd_outlet_node[i]];	
			a[j][node[sd_outlet_node[i]].n_n[l]] = a[j][sd_outlet_node[i]];
			e[j][node[sd_outlet_node[i]].n_n[l]] = e[j][sd_outlet_node[i]];
			//mu[j][node[sd_outlet_node[i]].n_n[l]] = mu[j][sd_outlet_node[i]];
			
			u[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = u[j][sd_outlet_node[i]];
			v[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = v[j][sd_outlet_node[i]];
			w[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = w[j][sd_outlet_node[i]];
			p[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = p[j][sd_outlet_node[i]];
			rho[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = rho[j][sd_outlet_node[i]];
			t[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = t[j][sd_outlet_node[i]];	
			a[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = a[j][sd_outlet_node[i]];
			e[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = e[j][sd_outlet_node[i]];
			//mu[j][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = mu[j][sd_outlet_node[i]];
	
			u[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][sd_outlet_node[i]];
			v[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][sd_outlet_node[i]];
			w[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][sd_outlet_node[i]];
			p[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][sd_outlet_node[i]];
			rho[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][sd_outlet_node[i]];
			t[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][sd_outlet_node[i]];	
			a[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][sd_outlet_node[i]];
			e[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][sd_outlet_node[i]];
			//mu[j][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_outlet_node[i]];
		}
		else if (node[sd_outlet_node[i]].loc == 7 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);

		}
		else if (node[sd_outlet_node[i]].loc == 8 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 9 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 10 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		
		else if (node[sd_outlet_node[i]].loc == 11 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 12 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 13 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 14 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);

		}
		else if (node[sd_outlet_node[i]].loc == 15 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 16 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 17 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 18 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		
		else if (node[sd_outlet_node[i]].loc == 19 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 20 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 21 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 22 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 23 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 24 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 25 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 26 )
		{
			u[j][sd_outlet_node[i]] = u[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			v[j][sd_outlet_node[i]] = v[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			w[j][sd_outlet_node[i]] = w[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			p[j][sd_outlet_node[i]] = p[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			t[j][sd_outlet_node[i]] = t[j][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			rho[j][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_outlet_node[i]]/t[j][sd_outlet_node[i]]);
			e[j][sd_outlet_node[i]] = p[j][sd_outlet_node[i]]/(0.4*rho[j][sd_outlet_node[i]]);
		}
	}

/****************************************************sd_inlet********************************************************/
	for (i=0; i< sd_inl_node; i++)
	{	
		if (node[sd_inlet_node[i]].loc > 0 && node[sd_inlet_node[i]].loc <= 6)
		{	
			if (node[sd_inlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_inlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_inlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_inlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_inlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_inlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}
		
			u[j][sd_inlet_node[i]]= 1.0;		
			v[j][sd_inlet_node[i]]= 0.0;	
			w[j][sd_inlet_node[i]]= 0.0;				
			p[j][sd_inlet_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[j][sd_inlet_node[i]]= 1.0;		
			rho[j][sd_inlet_node[i]]= (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
			a[j][sd_inlet_node[i]] = sqrt(1.4*p[j][sd_inlet_node[i]]/rho[j][sd_inlet_node[i]]); 

			u[j][node[sd_inlet_node[i]].n_n[l]] = u[j][sd_inlet_node[i]];
			v[j][node[sd_inlet_node[i]].n_n[l]] = v[j][sd_inlet_node[i]];
			w[j][node[sd_inlet_node[i]].n_n[l]] = w[j][sd_inlet_node[i]];
			p[j][node[sd_inlet_node[i]].n_n[l]] = p[j][sd_inlet_node[i]];
			rho[j][node[sd_inlet_node[i]].n_n[l]] = rho[j][sd_inlet_node[i]];
			t[j][node[sd_inlet_node[i]].n_n[l]] = t[j][sd_inlet_node[i]];	
			a[j][node[sd_inlet_node[i]].n_n[l]] = a[j][sd_inlet_node[i]];
			e[j][node[sd_inlet_node[i]].n_n[l]] = e[j][sd_inlet_node[i]];
			//mu[j][node[sd_inlet_node[i]].n_n[l]] = mu[j][sd_inlet_node[i]];
			
			u[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = u[j][sd_inlet_node[i]];
			v[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = v[j][sd_inlet_node[i]];
			w[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = w[j][sd_inlet_node[i]];
			p[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = p[j][sd_inlet_node[i]];
			rho[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = rho[j][sd_inlet_node[i]];
			t[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = t[j][sd_inlet_node[i]];	
			a[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = a[j][sd_inlet_node[i]];
			e[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = e[j][sd_inlet_node[i]];
			//mu[j][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = mu[j][sd_inlet_node[i]];
	
			u[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][sd_inlet_node[i]];
			v[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][sd_inlet_node[i]];
			w[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][sd_inlet_node[i]];
			p[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][sd_inlet_node[i]];
			rho[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][sd_inlet_node[i]];
			t[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][sd_inlet_node[i]];	
			a[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][sd_inlet_node[i]];
			e[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][sd_inlet_node[i]];
			//mu[j][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_inlet_node[i]];
		}
		else if (node[sd_inlet_node[i]].loc == 7 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 8 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 9 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 10 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		
		else if (node[sd_inlet_node[i]].loc == 11 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 12 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 13 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 14 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);

		}
		else if (node[sd_inlet_node[i]].loc == 15 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 16 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 17 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 18 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		
		else if (node[sd_inlet_node[i]].loc == 19 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 20 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 21 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 22 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 23 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 24 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 25 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 26 )
		{
			u[j][sd_inlet_node[i]] = u[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			v[j][sd_inlet_node[i]] = v[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			w[j][sd_inlet_node[i]] = w[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			p[j][sd_inlet_node[i]] = p[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			t[j][sd_inlet_node[i]] = t[j][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			rho[j][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[j][sd_inlet_node[i]]/t[j][sd_inlet_node[i]]);
			e[j][sd_inlet_node[i]] = p[j][sd_inlet_node[i]]/(0.4*rho[j][sd_inlet_node[i]]);
		}
	}
	
	/****************************************************sd_boundary********************************************************/
	for (i=0; i< sd_bou_node; i++)
	{	
		if (node[sd_boundary_node[i]].loc > 0 && node[sd_boundary_node[i]].loc <= 6)
		{	
			if (node[sd_boundary_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_boundary_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_boundary_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_boundary_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_boundary_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_boundary_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

/*			u[j][sd_boundary_node[i]]= 1.0;		
			v[j][sd_boundary_node[i]]= 0.0;	
			w[j][sd_boundary_node[i]]= 0.0;
			p[j][sd_boundary_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[j][sd_boundary_node[i]]= 1.0;		
			rho[j][sd_boundary_node[i]]= (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
			a[j][sd_boundary_node[i]] = sqrt(1.4*p[j][sd_boundary_node[i]]/rho[j][sd_boundary_node[i]]); 
*/
			u[j][sd_boundary_node[i]] = u[j][node[sd_boundary_node[i]].n_n[k]];
			v[j][sd_boundary_node[i]] = v[j][node[sd_boundary_node[i]].n_n[k]];
			w[j][sd_boundary_node[i]] = w[j][node[sd_boundary_node[i]].n_n[k]];
			p[j][sd_boundary_node[i]] = p[j][node[sd_boundary_node[i]].n_n[k]];
			rho[j][sd_boundary_node[i]] = rho[j][node[sd_boundary_node[i]].n_n[k]];
			t[j][sd_boundary_node[i]] = t[j][node[sd_boundary_node[i]].n_n[k]];	
			a[j][sd_boundary_node[i]] = a[j][node[sd_boundary_node[i]].n_n[k]];
			e[j][sd_boundary_node[i]] = e[j][node[sd_boundary_node[i]].n_n[k]];
			//mu[j][sd_boundary_node[i]] = mu[j][sd_boundary_node[i]];
			
			u[j][node[sd_boundary_node[i]].n_n[l]] = u[j][node[sd_boundary_node[i]].n_n[k]];
			v[j][node[sd_boundary_node[i]].n_n[l]] = v[j][node[sd_boundary_node[i]].n_n[k]];
			w[j][node[sd_boundary_node[i]].n_n[l]] = w[j][node[sd_boundary_node[i]].n_n[k]];
			p[j][node[sd_boundary_node[i]].n_n[l]] = p[j][node[sd_boundary_node[i]].n_n[k]];
			rho[j][node[sd_boundary_node[i]].n_n[l]] = rho[j][node[sd_boundary_node[i]].n_n[k]];
			t[j][node[sd_boundary_node[i]].n_n[l]] = t[j][node[sd_boundary_node[i]].n_n[k]];	
			a[j][node[sd_boundary_node[i]].n_n[l]] = a[j][node[sd_boundary_node[i]].n_n[k]];
			e[j][node[sd_boundary_node[i]].n_n[l]] = e[j][node[sd_boundary_node[i]].n_n[k]];
			//mu[j][node[sd_boundary_node[i]].n_n[l]] = mu[j][sd_boundary_node[i]];
			
			u[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = u[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			v[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = v[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			w[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = w[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			p[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = p[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			rho[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = rho[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			t[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = t[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];	
			a[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = a[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			e[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = e[j][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			//mu[j][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = mu[j][sd_boundary_node[i]];
	
			u[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];	
			a[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_boundary_node[i]];
			//mu[j][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_boundary_node[i]];
		}
		else if (node[sd_boundary_node[i]].loc == 7 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);

		}
		else if (node[sd_boundary_node[i]].loc == 8 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 9 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 10 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		
		else if (node[sd_boundary_node[i]].loc == 11 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 12 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 13 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 14 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);

		}
		else if (node[sd_boundary_node[i]].loc == 15 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 16 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 17 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 18 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		
		else if (node[sd_boundary_node[i]].loc == 19 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 20 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 21 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 22 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 23 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 24 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 25 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 26 )
		{
			u[j][sd_boundary_node[i]] = u[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			v[j][sd_boundary_node[i]] = v[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			w[j][sd_boundary_node[i]] = w[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			p[j][sd_boundary_node[i]] = p[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			t[j][sd_boundary_node[i]] = t[j][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			rho[j][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[j][sd_boundary_node[i]]/t[j][sd_boundary_node[i]]);
			e[j][sd_boundary_node[i]] = p[j][sd_boundary_node[i]]/(0.4*rho[j][sd_boundary_node[i]]);
		}		
	}
	
	/*****************************************************************************************************************************************************************/
	
	for (i=0; i< sd_wal_node; i++)
	{	
		if (node[sd_wall_node[i]].loc > 0 && node[sd_wall_node[i]].loc <= 6)
		{	
			if (node[sd_wall_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_wall_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_wall_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_wall_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_wall_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_wall_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[sd_wall_node[i]].n_n[k]];
			//p[j][sd_wall_node[i]]= 1.0/(1.4*Mach*Mach);
			t[j][sd_wall_node[i]]= t[j][node[sd_wall_node[i]].n_n[k]];	
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);		
			a[j][sd_wall_node[i]] = sqrt(1.4*p[j][sd_wall_node[i]]/rho[j][sd_wall_node[i]]);	
		
			u[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[j][node[sd_wall_node[i]].n_n[k]];
			v[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[j][node[sd_wall_node[i]].n_n[k]];
			w[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[j][node[sd_wall_node[i]].n_n[k]];
			p[j][node[sd_wall_node[i]].n_n[l]] = p[j][node[sd_wall_node[i]].n_n[k]];
			rho[j][node[sd_wall_node[i]].n_n[l]] = rho[j][node[sd_wall_node[i]].n_n[k]];
			t[j][node[sd_wall_node[i]].n_n[l]] = t[j][node[sd_wall_node[i]].n_n[k]];
			a[j][node[sd_wall_node[i]].n_n[l]] = a[j][node[sd_wall_node[i]].n_n[k]];
			e[j][node[sd_wall_node[i]].n_n[l]] = e[j][node[sd_wall_node[i]].n_n[k]];
			//mu[j][node[sd_wall_node[i]].n_n[l]] = mu[j][sd_wall_node[i]];
		
			u[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			v[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			w[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			p[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			rho[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			t[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			a[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			e[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[j][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			//mu[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[j][sd_wall_node[i]];
		
			u[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			a[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_wall_node[i]];
			
			
			/*******************************BLEED*********************************************/
		/* 	if ( (node[sd_wall_node[i]].x >= 1120.0 && node[sd_wall_node[i]].x <= 1220.0 && node[sd_wall_node[i]].y < 354.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1350.0 && node[sd_wall_node[i]].x <= 1410.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1015.238 && node[sd_wall_node[i]].x <= 1250.0 && node[sd_wall_node[i]].loc == 1) || (node[sd_wall_node[i]].x >= 1350.0 && node[sd_wall_node[i]].x <= 1482.0 && node[sd_wall_node[i]].loc == 1))
//			if ((node[sd_wall_node[i]].x >= 1226.5+700.0 && node[sd_wall_node[i]].x <= 1411.224+700.0 && node[sd_wall_node[i]].y < 674.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1411.224+700.0 && node[sd_wall_node[i]].x <= 1564.39+700.0 && node[sd_wall_node[i]].loc == 3))	
			{
			//	 if (node[sd_wall_node[i]].x >= 720.0 && node[sd_wall_node[i]].x <= 1226.5 && node[sd_wall_node[i]].y < 354.0 && node[sd_wall_node[i]].loc == 3)
			//	{
			//		Cd = Cdr1;
			//	} 
				 
				if (node[sd_wall_node[i]].x >= 1120.0 && node[sd_wall_node[i]].x <= 1220.0 && node[sd_wall_node[i]].y < 354.0 && node[sd_wall_node[i]].loc == 3)
				{
					Cd = Cdr2;
				}
				
				if (node[sd_wall_node[i]].x >= 1350.0 && node[sd_wall_node[i]].x <= 1410.0 && node[sd_wall_node[i]].loc == 3)
				{
					Cd = Cdr3;
				}
				
				if (node[sd_wall_node[i]].x >= 1015.238 && node[sd_wall_node[i]].x <= 1250.0 && node[sd_wall_node[i]].loc == 1)
				{
					Cd = Cdc1;
				}
		 	
				if (node[sd_wall_node[i]].x >= 1350.0 && node[sd_wall_node[i]].x <= 1482.0 && node[sd_wall_node[i]].loc == 1)
				{
					Cd = Cdc2;
				}
				
				if (node[sd_wall_node[i]].loc == 3)
				{
					Vn_sign = -1.0;
				}
				
				if (node[sd_wall_node[i]].loc == 1)
				{
					Vn_sign = 1.0;
				}
				
				rho[j][sd_wall_node[i]]= rho[j][node[sd_wall_node[i]].n_n[k]];
				angle = atan((node[node[sd_wall_node[i]].n_n[1]].y-node[sd_wall_node[i]].y)/(node[node[sd_wall_node[i]].n_n[1]].x-node[sd_wall_node[i]].x));
				angle2 = atan((node[node[node[sd_wall_node[i]].n_n[k]].n_n[1]].y-node[node[sd_wall_node[i]].n_n[k]].y)/(node[node[node[sd_wall_node[i]].n_n[k]].n_n[1]].x-node[node[sd_wall_node[i]].n_n[k]].x));
				determinant = cos(angle)*cos(angle)+sin(angle)*sin(angle);

				Vt = u[j][node[sd_wall_node[i]].n_n[k]]*cos(angle)+v[j][node[sd_wall_node[i]].n_n[k]]*sin(angle);
				Vn = Vn_sign*Cd*sqrt(1.4*p[j][node[sd_wall_node[i]].n_n[k]]/rho[j][node[sd_wall_node[i]].n_n[k]]);
				Vn2 = v[j][node[sd_wall_node[i]].n_n[k]]*cos(angle2)-u[j][node[sd_wall_node[i]].n_n[k]]*sin(angle2);
				
				p[j][sd_wall_node[i]] = p[j][node[sd_wall_node[i]].n_n[k]]-2.0*rho[j][sd_wall_node[i]]*Vn*(Vn-Vn2);				
				u[j][sd_wall_node[i]] = (cos(angle)/determinant)*Vt-(sin(angle)/determinant)*Vn;
				v[j][sd_wall_node[i]] = (sin(angle)/determinant)*Vt+(cos(angle)/determinant)*Vn;
				w[j][sd_wall_node[i]]= 0.0;
				t[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/rho[j][sd_wall_node[i]]);
				e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);		
				a[j][sd_wall_node[i]] = sqrt(1.4*p[j][sd_wall_node[i]]/rho[j][sd_wall_node[i]]);	
				
				u[j][node[sd_wall_node[i]].n_n[l]] = u[j][sd_wall_node[i]];
				v[j][node[sd_wall_node[i]].n_n[l]] = v[j][sd_wall_node[i]];
				w[j][node[sd_wall_node[i]].n_n[l]] = w[j][sd_wall_node[i]];
				p[j][node[sd_wall_node[i]].n_n[l]] = p[j][sd_wall_node[i]];
				rho[j][node[sd_wall_node[i]].n_n[l]] = rho[j][sd_wall_node[i]];
				t[j][node[sd_wall_node[i]].n_n[l]] = t[j][sd_wall_node[i]];	
				a[j][node[sd_wall_node[i]].n_n[l]] = a[j][sd_wall_node[i]];
				e[j][node[sd_wall_node[i]].n_n[l]] = e[j][sd_wall_node[i]];
				
				u[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = u[j][sd_wall_node[i]];
				v[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = v[j][sd_wall_node[i]];
				w[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = w[j][sd_wall_node[i]];
				p[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[j][sd_wall_node[i]];
				rho[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[j][sd_wall_node[i]];
				t[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[j][sd_wall_node[i]];	
				a[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[j][sd_wall_node[i]];
				e[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[j][sd_wall_node[i]];
		
				u[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[j][sd_wall_node[i]];
				v[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[j][sd_wall_node[i]];
				w[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[j][sd_wall_node[i]];
				p[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][sd_wall_node[i]];
				rho[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][sd_wall_node[i]];
				t[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][sd_wall_node[i]];	
				a[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][sd_wall_node[i]];
				e[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][sd_wall_node[i]];
				
			} 		
			 */ 	
			
			
		}
		else if (node[sd_wall_node[i]].loc == 7 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);

		}
		else if (node[sd_wall_node[i]].loc == 8 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 9 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 10 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		
		else if (node[sd_wall_node[i]].loc == 11 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 12 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 13 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 14 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 15 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[sd_wall_node[i]].n_n[5]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[5]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 16 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[sd_wall_node[i]].n_n[1]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[1]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 17 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]]= p[j][node[node[sd_wall_node[i]].n_n[4]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[4]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 18 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[3]].n_n[0]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[3]].n_n[0]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		
		else if (node[sd_wall_node[i]].loc == 19 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[5]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[5]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 20 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[1]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[1]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 21 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[4]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[4]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 22 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[3]].n_n[2]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[3]].n_n[2]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 23 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[1]].n_n[5]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[1]].n_n[5]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 24 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[1]].n_n[4]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[1]].n_n[4]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 25 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 26 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 29 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 30 )
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			t[j][sd_wall_node[i]] = t[j][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);
		}
		
		
		if (node[sd_wall_node[i]].ID == 10 && node[sd_wall_node[i]].loc != 29 && node[sd_wall_node[i]].loc != 30)
		{
			u[j][sd_wall_node[i]]= 0.0;
			v[j][sd_wall_node[i]]= 0.0;
			w[j][sd_wall_node[i]]= 0.0;
			p[j][sd_wall_node[i]] = p[j][node[sd_wall_node[i]].n_n[3]];
			t[j][sd_wall_node[i]] = t[j][node[sd_wall_node[i]].n_n[3]];
			rho[j][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[j][sd_wall_node[i]]/t[j][sd_wall_node[i]]);
			e[j][sd_wall_node[i]] = p[j][sd_wall_node[i]]/(0.4*rho[j][sd_wall_node[i]]);		
			
			for (h=0;h<3;h++)
			{
				if (h==0)
				{
					l = 2; /********no element on SOUTH**************/
					k = 0;
				}
				if (h==1) 
				{
					l = 0; /********no element on NORTH**************/
					k = 2;
				}
				if (h==2)
				{
					l = 1; /********no element on EAST**************/
					k = 3;										
				}
		
				if (h<=1)
				{
					u[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[j][singular[sd_wall_node[i]].n_n[k]];
					v[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[j][singular[sd_wall_node[i]].n_n[k]];
					w[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[j][singular[sd_wall_node[i]].n_n[k]];
					p[j][node[sd_wall_node[i]].n_n[l]] = p[j][singular[sd_wall_node[i]].n_n[k]];
					rho[j][node[sd_wall_node[i]].n_n[l]] = rho[j][singular[sd_wall_node[i]].n_n[k]];
					t[j][node[sd_wall_node[i]].n_n[l]] = t[j][singular[sd_wall_node[i]].n_n[k]];	
					a[j][node[sd_wall_node[i]].n_n[l]] = a[j][sd_wall_node[i]];
					e[j][node[sd_wall_node[i]].n_n[l]] = e[j][singular[sd_wall_node[i]].n_n[k]];
					//mu[j][node[sd_wall_node[i]].n_n[l]] = mu[j][sd_wall_node[i]];
							
					u[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					v[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					w[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					p[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					rho[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					t[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];	
					a[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[j][sd_wall_node[i]];
					e[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					//mu[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[j][sd_wall_node[i]];
		
					u[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][sd_wall_node[i]];
					e[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_wall_node[i]];	
				}
				if (h == 2)
				{
					u[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[j][singular[sd_wall_node[i]].n_n[k]];
					v[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[j][singular[sd_wall_node[i]].n_n[k]];
					w[j][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[j][singular[sd_wall_node[i]].n_n[k]];
					p[j][node[sd_wall_node[i]].n_n[l]] = p[j][singular[sd_wall_node[i]].n_n[k]];
					rho[j][node[sd_wall_node[i]].n_n[l]] = rho[j][singular[sd_wall_node[i]].n_n[k]];
					t[j][node[sd_wall_node[i]].n_n[l]] = t[j][singular[sd_wall_node[i]].n_n[k]];	
					a[j][node[sd_wall_node[i]].n_n[l]] = a[j][sd_wall_node[i]];
					e[j][node[sd_wall_node[i]].n_n[l]] = e[j][singular[sd_wall_node[i]].n_n[k]];
					//mu[j][node[sd_wall_node[i]].n_n[l]] = mu[j][sd_wall_node[i]];
				
					u[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					v[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					w[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					p[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					rho[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					t[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];	
					a[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[j][sd_wall_node[i]];
					e[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[j][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					//mu[j][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[j][sd_wall_node[i]];
		
					u[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[j][sd_wall_node[i]];
					e[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[j][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[j][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[j][sd_wall_node[i]];	
				}
			}
		}
	}
	
}







