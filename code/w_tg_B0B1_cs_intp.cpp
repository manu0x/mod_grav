using namespace std;

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <gsl/gsl_sf_bessel.h>

#include "./other_source/spline.c"
#include "./other_source/power.c"


#define c_box 2.99
#define twopie 2.0*M_PI
#define H0byc (0.001/2.99)




class cosmo_bigravity
{

	private:
	double  ratio,b,beta,gamma;
	int model;

	public:
	double omega_dm_0, H0,B1;
	double D1_0,net_chi; 
	double  B0;
	double dai;
	

	double *cs_D1,*cs_f,*cs_cmd,*a_sp,*D1_sp,*f_sp,*cmd_sp;
	cosmo_bigravity(double omega_dm_0_val=0.3,double B1param=0.0,double h_val=0.7,int model=2)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = (h_val/c_box)*0.001;
		model = model;
		B1 = B1param;
		B0 = 3.0*(1.0-omega_dm_0-B1*B1/3.0);
		ratio = (1.0-omega_dm_0)/(omega_dm_0);
		printf("H0 %lf  B1  %lf\n",H0,B1);

		

		TFmdm_set_cosm(omega_dm_0, 0.0223/(h_val*h_val), 0.0,
		1, (1.0-omega_dm_0), h_val, 1.0);

	}

	double Hsqr(double a,double a0=1.0)
	{
		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		val = H0*H0*(   0.5*omega + B0/6.0  + sqrt((0.5*omega + B0/6.0)*(0.5*omega + B0/6.0) + B1*B1/3.0 ) );
		return(val);


	}


	void g(double x,double dr[4])
	{
 	   double b0t,b1t;
 	   b0t = B0/6.0 + 0.5*x;
    	   b1t = B1*B1/3.0;
    	   dr[0] =  b0t + sqrt(b0t*b0t + b1t);
    	   dr[1] = 0.5 + 0.5*b0t/sqrt(b1t + b0t*b0t);
    	   dr[2] = -0.25*b0t*b0t/pow(b1t + b0t*b0t,1.5) + 0.25/sqrt(b1t + b0t*b0t);
    	   dr[3] = (3.0*b0t*b0t*b0t/pow(b1t + b0t*b0t,2.5) -3.0*b0t/pow(b1t + b0t*b0t,1.5))/8.0;	



	}


	double Hb_t_by_Hb2(double a,double a0=1.0)
	{

		double val,omega_dm;
		omega_dm = omega_dm_0*pow(a0/a,3.0);		
		beta = 0.5*omega_dm + B0/6.0;
		
		val = -(3.0/4.0)*omega_dm/sqrt(beta*beta + B1*B1/3.0);
		return(val);



	}

	double H_Diff(double a, double delta, double a0=1.0)
	{
		double diff, term1,term2,omega_dm;
	
		omega_dm = omega_dm_0*pow(a0/a,3.0);		
		beta = 0.5*omega_dm + B0/6.0;
		b = 0.5*omega_dm*(1.0+delta) + B0/6.0;
		gamma = (b+ sqrt(b*b + B1*B1/3.0))/(beta+ sqrt(beta*beta + B1*B1/3.0));
		term1 = 1.0/sqrt(beta*beta + B1*B1/3.0);
		term2 = (1.0+delta)*gamma/sqrt(b*b + B1*B1/3.0);
		diff = (gamma-1.0) + (3.0/4.0)*omega_dm*(term1-term2) ;

		//printf("diff %lf  bigr  %lf\n",diff,1.0+delta/(1.0+ratio*pow(a/a0,3.0)));
		return(diff);
		

	}

	double delta_aa(double a, double delta, double delta_a)
	{
		double HbtbyHb2, diff, acc;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		diff = H_Diff(a, delta);
		
		
		acc = -(1.0/a)*(3.0 + HbtbyHb2)*delta_a + (4.0/3.0)*delta_a*delta_a/(1.0+delta) - 3.0*(1.0+delta)*diff/(a*a);
		
	
		return(acc);

	}


	double lin_delta_aa(double a, double delta, double delta_a,double a0=1.0)
	{
		double HbtbyHb2, diff, acc, theta,kappa,lambda;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		diff = H_Diff(a, delta);
		
		theta = pow(a0/a,3.0);
		kappa = -3.0*( 9.0*omega_dm_0*theta + 6.0*B1*B1*omega_dm_0*theta + pow(B1,4.0)*theta*omega_dm_0 
				-18.0*theta*omega_dm_0*omega_dm_0 + 6.0*B1*B1*theta*omega_dm_0*omega_dm_0 - 9.0*theta*theta*omega_dm_0*omega_dm_0 
				+ 3.0*B1*B1*theta*theta*omega_dm_0*omega_dm_0 + 9.0*theta*omega_dm_0*omega_dm_0*omega_dm_0
				+ 9.0*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0 - 18.0*omega_dm_0*omega_dm_0*omega_dm_0*theta*theta*theta
				+ 54.0*theta*theta*omega_dm_0*omega_dm_0*sqrt(B1*B1/3.0 + (B0/6.0 + 0.5*omega_dm_0*theta)*(B0/6.0 + 0.5*omega_dm_0*theta)));
		lambda = (9.0 + 6.0*B1*B1 + B1*B1*B1*B1 -18.0*omega_dm_0 + 6.0*B1*B1*omega_dm_0 + 18.0*theta*omega_dm_0 
				-6.0*B1*B1*theta*omega_dm_0 + 9.0*omega_dm_0*omega_dm_0 - 18.0*theta*omega_dm_0*omega_dm_0
				+ 9.0*theta*theta*omega_dm_0*omega_dm_0);
		lambda = 2.0*pow(lambda,1.5);
		
		diff = kappa*delta/lambda;
		acc = -(1.0/a)*(3.0 + HbtbyHb2)*delta_a  - 3.0*diff/(a*a);
	
		return(acc);

	}

	double delta_lin_lue(double a, double delta, double delta_a,int de_mode,double a0=1.0)
	{
		double x,acc,gv[4];

		x =  omega_dm_0*pow(a0/a,3.0);
		  
		  g(x,gv);
		 
		if(de_mode)
		{
			acc = -(1.0-0.5*x*(gv[1]/gv[0]))*3.0*(delta_a/a) + 1.5*x*delta/(gv[0]*a*a);

		}

		else
		{
			acc = -(1.0-0.5*x*(gv[1]/gv[0]))*3.0*(delta_a/a) + 1.5*x*(gv[1]+ 3.0*x*gv[2])*delta/(gv[0]*a*a);

		}

		
	
		return(acc);

	}



	void pert_delta_aa(double *acc, double D[5], double D_a[5],double a,double a0=1.0)
	{
		double HbtbyHb2, diff, theta,kappa1,lambda1,kappa2,lambda2;	
		double c1,c2;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		
		
		theta = pow(a0/a,3.0);
		kappa1 = -3.0*( 9.0*omega_dm_0*theta + 6.0*B1*B1*omega_dm_0*theta + pow(B1,4.0)*theta*omega_dm_0 
				-18.0*theta*omega_dm_0*omega_dm_0 + 6.0*B1*B1*theta*omega_dm_0*omega_dm_0 - 9.0*theta*theta*omega_dm_0*omega_dm_0 
				+ 3.0*B1*B1*theta*theta*omega_dm_0*omega_dm_0 + 9.0*theta*omega_dm_0*omega_dm_0*omega_dm_0
				+ 9.0*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0 - 18.0*omega_dm_0*omega_dm_0*omega_dm_0*theta*theta*theta
				+ 54.0*theta*theta*omega_dm_0*omega_dm_0*sqrt(B1*B1/3.0 + (B0/6.0 + 0.5*omega_dm_0*theta)*(B0/6.0 + 0.5*omega_dm_0*theta)));
		lambda1 = (9.0 + 6.0*B1*B1 + B1*B1*B1*B1 -18.0*omega_dm_0 + 6.0*B1*B1*omega_dm_0 + 18.0*theta*omega_dm_0 
				-6.0*B1*B1*theta*omega_dm_0 + 9.0*omega_dm_0*omega_dm_0 - 18.0*theta*omega_dm_0*omega_dm_0
				+ 9.0*theta*theta*omega_dm_0*omega_dm_0);
		lambda1 = 2.0*pow(lambda1,1.5);
		
		c1 = kappa1/lambda1;

		kappa2 = 27.0*( -36.0*B1*B1*theta*theta*omega_dm_0*omega_dm_0 -24.0*B1*B1*B1*B1*theta*theta*omega_dm_0*omega_dm_0 
				-4.0*pow(B1,6.0)*theta*theta*omega_dm_0*omega_dm_0 + 72.0*B1*B1*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0
				-24.0*pow(B1,4.0)*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0 + 9.0*B1*B1*theta*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0
				- 3.0*pow(B1,4.0)*theta*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0
				-36.0*B1*B1*theta*theta*omega_dm_0*omega_dm_0*omega_dm_0*omega_dm_0 - 9.0*B1*B1*theta*theta*theta*pow(omega_dm_0,4.0)
				+ 45.0*B1*B1*pow(theta,4.0)*pow(omega_dm_0,4.0));
		lambda2 = 9.0+6.0*B1*B1+pow(B1,4.0) - 18.0*omega_dm_0 + 6.0*B1*B1*omega_dm_0 
				+ 18.0*theta*omega_dm_0 - 6.0*B1*B1*theta*omega_dm_0 + 9.0*omega_dm_0*omega_dm_0
				- 18.0*theta*omega_dm_0*omega_dm_0 + 9.0*theta*theta*omega_dm_0*omega_dm_0;
		lambda2 = pow(lambda2,2.5);
		lambda2 = lambda2*(3.0 - B1*B1 - 3.0*omega_dm_0 + 3.0*omega_dm_0*theta + 6.0*sqrt(B1*B1/3.0 +
						 (B0/6.0 + 0.5*omega_dm_0*theta)*(B0/6.0 + 0.5*omega_dm_0*theta) ));
	
		c2 = kappa2/lambda2;

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);
		acc[2] = delta_aa(a, D[2],  D_a[2]);
		acc[3] = delta_lin_lue(a,D[3],D_a[3],1);
		acc[4] = delta_lin_lue(a,D[4],D_a[4],0);




	}



     int run_cosmo(FILE *fp,FILE *fp_de,FILE *fp_mg,int de_mode=1,double da=0.000001)
     {

	 double D[5],D_a[5],D_rk[5][5], D_a_rk[5][5], acc[5];
	 double D_i[5], D_a_i[5];

	 double a,ai,ai_burn,a0,ak;
	 double rk_coef[4] = {0.5,0.5,1.0,1.0};
	 
	 double isw_potn,isw_potn_de,isw_potn_mg,fv;
	
	double in_z,chi1,chi2,chi3,chi,cur_z,gv[4],x;
	double de_intzc,de_intzc1;
	double mg_intzc,mg_intzc1;
	

	int i,j,theory,acntr,aN;

	

	


	dai = da;
	ai = 0.000908;
	ai_burn = ai*0.1;
	a0 = 1.0;

	aN = (int)((a0-ai)/da) + 1;


	D_a_i[0] = 1.0;
	D_i[0] = ai_burn*D_a_i[0];
	
	D[0] = D_i[0];
	D_a[0] = D_a_i[0];

	D_i[1] = (34.0/7.0)*D_i[0]*D_i[0]/3.0;
	D_a_i[1] = 2.0*(34.0/7.0)*D_i[0]*D_a_i[0]/3.0;//D_i[1]/ai;
	
	D[1] = D_i[1];
	D_a[1] = D_a_i[1];

	D_i[2] = D_i[0];
	D_a_i[2] = D_a_i[0];

	D_i[3] = D_i[0];
	D_a_i[3] = D_a_i[0];

	D_i[4] = D_i[0];
	D_a_i[4] = D_a_i[0];
	
	D[2] = D_i[2];
	D_a[2] = D_a_i[2];


	D[3] = D_i[3];
	D_a[3] = D_a_i[3];

	D[4] = D_i[4];
	D_a[4] = D_a_i[4];


	a_sp = new double[aN];
	D1_sp = new double[aN];
	f_sp = new double[aN];
	cmd_sp = new double[aN];

	int burn = 1,cntr, cntr_spl;


	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);

	for(a=ai_burn,cntr=0;a<=a0;a+=da,++cntr)
	{
		ak = a;


		

		if((a>ai)&&(burn))
		{
			burn = 0;
			acntr = 0;
			cntr_spl=0;

			D_i[0] = D[0];
			D_a_i[0] = D_a[0];
	
			D_i[2] = D[2];
			D_a_i[2] = D_a[2];

			//D[1] = 4.858*D[0]*D[0]/3.0;
			//D_a[1] = 4.858*2.0*D_a[0]*D[0]/3.0;
		
			//D_i[1] = D[1];
			//D_a_i[1] = D_a[1];
		
			//printf("intial dc is %lf  %lf  %lf\n",D_i[0],D_i[1],D_i[2]);
			
			chi = 0.0;
			cntr_spl=0;


		}



		D_rk[0][0] = D[0];
		D_a_rk[0][0] = D_a[0]; 

		D_rk[0][1] = D[1];
		D_a_rk[0][1] = D_a[1]; 

		D_rk[0][2] = D[2];
		D_a_rk[0][2] = D_a[2]; 

		D_rk[0][3] = D[3];
		D_a_rk[0][3] = D_a[3]; 

		D_rk[0][4] = D[4];
		D_a_rk[0][4] = D_a[4]; 
	
		
		//
		//if(((cntr%1000)==0))
		if(!(burn)&&((cntr%1000)==0))
		{
			x =  omega_dm_0*pow(a0/a,3.0);
		
		        g(x,gv);	

			fv = a*(D_a[0]/D[0]);		
			isw_potn = D[0]*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) );

			fv = a*(D_a[3]/D[3]);
			isw_potn_de = (1.0-fv)*D[3];

			fv = a*(D_a[4]/D[4]);
			isw_potn_mg = D[4]*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) );
	
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				a,D[0]/a,D[1],D[2],D[0]*ai/(a*D_i[0]),D[1]*ai/(a*D_i[1]),D[0]/D_i[0],D[1]/D_i[1],D[2]/D_i[2],
				3.0*D[1]/(D[0]*D[0]),3.0*D[2]/(D[0]*D[0]),gv[0],isw_potn);

			fprintf(fp_de,"%lf\t%lf\t%lf\t%lf\t\n",
				a,D[3]/a,isw_potn_de,de_intzc);

			fprintf(fp_mg,"%lf\t%lf\t%lf\t%lf\t\n",
				a,D[4]/a,isw_potn_mg,mg_intzc);



		}
	

		if(a>=ai)
		{
		 	cmd_sp[cntr_spl] = chi;
			D1_sp[cntr_spl] = D[0];
			f_sp[cntr_spl] =  a*(D_a[0]/D[0]);
			a_sp[cntr_spl] = a;
			cntr_spl++;
		
	
		  x =  omega_dm_0*pow(a0/a,3.0);
		  chi1 = 1.0/(a*a*sqrt(Hsqr(a)));
		  chi2 = 1.0/((a+0.5*da)*(a+0.5*da)*sqrt(Hsqr(a+0.5*da)));
		  chi3 = 1.0/((a+da)*(a+da)*sqrt(Hsqr(a+da)));
		  chi+= (da/6.0)*(chi1+4.0*chi2+chi3)*(H0/H0byc);
		  g(x,gv);


		  fv = a*(D_a[3]/D[3]);
		  de_intzc1 = sqrt(gv[0])*(D[3]/a)*(D[3]/a)*( (1.0-fv));

		  fv = a*(D_a[4]/D[4]);
		  mg_intzc1 = sqrt(gv[0])*(D[4]/a)*(D[4]/a)*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) );

		  if(cntr==0)
		  {
   			
			de_intzc+= de_intzc1*da/3.0;	
			mg_intzc+= mg_intzc1*da/3.0;
		  }

		 else
		 {
			if(cntr%2==0)
			{
				de_intzc+= 2.0*de_intzc1*da/3.0;	
				mg_intzc+= 2.0*mg_intzc1*da/3.0;
			
			}
			

			else
			{
				de_intzc+= 4.0*de_intzc1*da/3.0;	
				mg_intzc+= 4.0*mg_intzc1*da/3.0;
			
			}


		 }

		 
		 
		 // cur_z = sqrt(gv[0])*D[0]*D[0]*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) )*wgv*j0;

		}
		

		for(i=1;i<=4;++i)
		{
			
			pert_delta_aa(acc, D_rk[0], D_a_rk[0],ak);

			
		   for(j=0;j<5;++j)
			{D_a_rk[i][j] = da*acc[j];
			 D_rk[i][j] = da*D_a_rk[0][j];
				

			 

			 D_rk[0][j] = D[j] + rk_coef[i-1]*D_rk[i][j];
			 D_a_rk[0][j] = D_a[j] + rk_coef[i-1]*D_a_rk[i][j];
			}

			ak = a + rk_coef[i-1]*da;
			if(a==ai)
			 printf("acc %.10lf  %.10lf\n",D_rk[i][0],D_rk[i][1]);

		}


	    for(j=0;j<5;++j)
		{
		 D[j] = D[j] + (1.0/6.0)*(D_rk[1][j]+2.0*D_rk[2][j]+2.0*D_rk[3][j]+D_rk[4][j]);
		 D_a[j] = D_a[j] + (1.0/6.0)*(D_a_rk[1][j]+2.0*D_a_rk[2][j]+2.0*D_a_rk[3][j]+D_a_rk[4][j]);



		

		}


	   }

	x =  omega_dm_0*pow(a0/a,3.0);
		
	g(x,gv);	
	fv = a*(D_a[0]/D[0]);	

	printf("gv0  %lf\n",gv[0]);	

	isw_potn = D[0]*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) );

	fv = a*(D_a[3]/D[3]);
	isw_potn_de = (1.0-fv)*D[3];

	fv = a*(D_a[4]/D[4]);
	isw_potn_mg = D[4]*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) );
	
	fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				a,D[0],D[1],D[2],D[0]*ai/(a*D_i[0]),D[1]*ai/(a*D_i[1]),D[0]/D_i[0],D[1]/D_i[1],D[2]/D_i[2],
				3.0*D[1]/(D[0]*D[0]),3.0*D[2]/(D[0]*D[0]),gv[0],isw_potn);

	fprintf(fp_de,"%lf\t%lf\t%lf\t%lf\t\n",
				a,D[3]/a,isw_potn_de,de_intzc);

	fprintf(fp_mg,"%lf\t%lf\t%lf\t%lf\t\n",
				a,D[4]/a,isw_potn_mg,mg_intzc);
	D1_0  = D[0];
	net_chi = chi;

	cs_D1 = new double[3*aN];
	cs_f = new double[3*aN];
	cs_cmd = new double[3*aN];

	printf("ffff %d\n",aN);
	double cs_D1s[aN][3],cs_fs[aN][3],cs_cmds[aN][3];

	

	int chk1,chk2,chk3;
	
	chk1 = spline(a_sp,D1_sp,aN,cs_D1s);
	chk2 = spline(a_sp,f_sp,aN,cs_fs);
	chk3 = spline(a_sp,cmd_sp,aN,cs_cmds);


	for(i=0;i<aN;++i)
	{

		cs_D1[3*i] = cs_D1s[i][0];
		cs_f[3*i] = cs_fs[i][0];
		cs_cmd[3*i] = cs_cmds[i][0];


		cs_D1[3*i+1] = cs_D1s[i][1];
		cs_f[3*i+1] = cs_fs[i][1];
		cs_cmd[3*i+1] = cs_cmds[i][1];

		cs_D1[3*i+2] = cs_D1s[i][2];
		cs_f[3*i+2] = cs_fs[i][2];
		cs_cmd[3*i+2] = cs_cmds[i][2];


	}

	if((!chk1)&&(!chk2)&&(!chk3))
	printf("Spline no err and D1_0 is %lf\n",D1_0);
	else
	printf("Errrr in spline %d %d %d");



	if(isnan(D[0]+D[1]+D[2]))
	return 0;
	else
	return 1;

	}

	
	

};




double prim_R(double k,double h=0.67)
{
	double pr;
	double As,kp,ns;
	As = 2.101e-9;
	kp= 0.05/h;
	ns=1.0;

	pr = (2.0*M_PI*M_PI/(k*k*k))*As*pow(k/kp,ns-1.0);

	return(pr);



}

double TF(double k)
{

	double tfv;
	tfv = TFmdm_onek_hmpc(k);
	return(tfv);
}



double wg(double z)
{

	double n0,z0,wgv;
	z0 = 0.49/3.0;
	n0 = (1.18)*(1e5);
	wgv = 0.5*z*z*exp(-z/z0)/pow(z0,3.0);
	return(wgv);
}

double Pk(double k,double om_0=0.3,double D1=1.0)
{
	double prv,tfv,pv;
	prv = prim_R(k);
	tfv = TF(k);
	pv = (4.0/25.0)*(k*k/om_0)*(k*k/om_0)*prv*tfv*tfv*D1*D1/pow(H0byc,4.0);

	return (pv);

}




void test_pk(double omega_dm_0,double D1_0)
{
	double pkv,primkv,tfkv,k,lk;
	FILE *fptest= fopen("testpk2.txt","w");
	for(lk=-12.0;lk<=1.0;lk+=0.05)
	{
		k=exp(lk);
		pkv = Pk(k,omega_dm_0,D1_0);
		primkv = prim_R(k);
		tfkv = TF(k);
		fprintf(fptest,"%lf\t%lf\t%lf\t%lf\n",k,pkv,primkv,tfkv);


	}


}

double int_z(int argc,char *argv[],double theta,double k,cosmo_bigravity bimet,double da =  0.0001)
{


  int i,j,cur_i;
  double in_z,D1,chi,gv[4],fv,j0,wgv,ktc,x;
  double ai,a0,a,delta;
  double cur_int;
  double as[3];
  
 
  double dsmp[3]={0.0,0.5*da,da};
  double smp[3]={1.0,4.0,1.0};

  ai = 0.001;
  a0 = 1.0;

  in_z = 0.0;
  
  for(a=ai;a<=a0;a+=da)
  {

	

	for(j=0;j<3;++j)
	{
		
		as[j] = a+dsmp[j];
		cur_i = (int)((as[j]-bimet.a_sp[0])/bimet.dai);
		delta = (a-bimet.a_sp[cur_i]);
		chi = bimet.cmd_sp[cur_i]+bimet.cs_cmd[3*cur_i]*delta+bimet.cs_cmd[3*cur_i+1]*delta*delta+bimet.cs_cmd[3*cur_i+2]*delta*delta*delta;	
		fv = bimet.f_sp[cur_i]+bimet.cs_f[3*cur_i]*delta+bimet.cs_f[3*cur_i+1]*delta*delta+bimet.cs_f[3*cur_i+2]*delta*delta*delta;
		D1 = bimet.D1_sp[cur_i]+bimet.cs_D1[3*cur_i]*delta+bimet.cs_D1[3*cur_i+1]*delta*delta+bimet.cs_D1[3*cur_i+2]*delta*delta*delta;

		x =  bimet.omega_dm_0*pow(a0/as[j],3.0);
		bimet.g(x,gv);
		
		ktc = k*theta*(bimet.net_chi-chi);
		j0 = gsl_sf_bessel_J0(ktc);
		wgv = wg(a0/a-1.0);

		cur_int = sqrt(gv[0])*(D1/as[j])*(D1/as[j])*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) )*wgv*j0;
		//cur_int = sqrt(gv[0])*(D1/as[j])*(D1/as[j])*( (1.0-fv))*wgv*j0;

		in_z+=(da/6.0)*smp[j]*cur_int;
			
		

	}

	


  }

  

  return(in_z);



}






double int_k(int argc,char *argv[],double theta,cosmo_bigravity bimet,double kini=0.0001,double kend=1.0,double dk=0.01)
{

	double in_k,in_k1,in_k2,in_k3;
	double pk1,pk2,pk3;
	double k;
	in_k = 0.0;

	double D1_0 = bimet.D1_0;
	
	printf("\n\nTHETA %lf\n",theta);
	
	for(k=kini;k<=kend;k+=dk)
	{  

	   pk1 = Pk(k,bimet.omega_dm_0,D1_0);
	   pk2 = Pk(k+0.5*dk,bimet.omega_dm_0,D1_0);
	   pk3 = Pk(k+dk,bimet.omega_dm_0,D1_0);

	   in_k1 = int_z(argc,argv,theta,k,bimet);//printf("theta %lf  k_frac %lf  %lf\n",theta,k,in_k1);
	   in_k1 = pk1*in_k1/k;

	   
	
	   in_k2 = int_z(argc,argv,theta,k+0.5*dk,bimet);
	   in_k2 = pk2*in_k2/(k+0.5*dk);

	   in_k3 = int_z(argc,argv,theta,k+dk,bimet);
	   in_k3 = pk3*in_k3/(k+dk);

	   in_k+= (dk/6.0)*(in_k1+4.0*in_k2+in_k3);


	}
	
	return(in_k);


}


double int_logk(int argc,char *argv[],double theta,cosmo_bigravity bimet,double lkini=-12.0,double lkend=1.0,double dlk=0.08)
{

	double in_k,in_k1,in_k2,in_k3;
	double pk1,pk2,pk3;
	double k,lk,k1,k2,k3;
	in_k = 0.0;

	double D1_0 = bimet.D1_0;
	
	//printf("\n\nTHETA %lf\n",theta);
	
	for(lk=lkini;lk<=lkend;lk+=dlk)
	{  
	
	   k1 = exp(lk);
	   k2 = exp(lk+0.5*dlk);
	   k3 = exp(lk+dlk);

	   pk1 = Pk(k1,bimet.omega_dm_0,D1_0);
	   pk2 = Pk(k2,bimet.omega_dm_0,D1_0);
	   pk3 = Pk(k3,bimet.omega_dm_0,D1_0);


	   in_k1 = int_z(argc,argv,theta,k1,bimet);//printf("theta %lf  k_frac %lf  %lf\n",theta,k,in_k1);
	   in_k1 = pk1*in_k1;

	   
	
	   in_k2 = int_z(argc,argv,theta,k2,bimet);
	   in_k2 = pk2*in_k2;

	   in_k3 = int_z(argc,argv,theta,k3,bimet);
	   in_k3 = pk3*in_k3;

	   in_k+= (dlk/6.0)*(in_k1+4.0*in_k2+in_k3);


	}
	
	return(in_k);


}





//int f_sigma_cal(int argc,char *argv[],double sig8,double *al,double *fs8,int n,int ax=1)



int main(int argc,char *argv[])
{
	double sig8  = 0.79;
	double a_l[12] = {0.02,0.08,0.15,0.19,0.23,0.35,0.45,0.55,0.65,0.78,0.89,0.999};
	double z_l[24] = {0.001,0.08,0.15,0.19,0.23,0.35,0.45,0.55,0.65,0.78,0.89,0.999,1.001,1.08,1.15,1.19,1.23,1.35,1.45,1.55,1.65,1.78,1.89,1.999};
	int n = 24,b;
	int ax = 0;
	double *data;


	double om_dm_0, model_param;
	double D1_0; int theory;

	double multi_fac,T0,bias;
	double theta,thetai,thetaend,dtheta,theta_r,wgt;

	
	thetai = 0.05;
	thetaend = 150.0;
	dtheta = 0.4*exp(1.0);

	string fname = "delta";
	string fname2 = "DE";
	string fname3 = "MG";
	string fname4 = "wgt";
	string argstr = argv[1];
	string om_str = argv[2];
	string mod_str = argv[3];
	string extstr = ".txt";
	string us = "_";

	
	
	om_dm_0 = atof(argv[2]);
	model_param = atof(argv[3]);

	T0 = 2.72548;
	bias = 5.47;

	multi_fac = 3.0*T0*twopie*twopie*bias*om_dm_0*H0byc*H0byc*H0byc;
	//multi_fac = 3.0*T0*bias*om_dm_0*H0byc*H0byc*H0byc/twopie;

	

	cosmo_bigravity cosmo_model_bigravity(om_dm_0,model_param);

	
	

	

	fname = fname+us+argstr+us+"om"+us+om_str+us+"B1"+us+mod_str;
	fname2 = fname2+us+argstr+us+"om"+us+om_str+us+"B1"+us+mod_str;
	fname3 = fname3+us+argstr+us+"om"+us+om_str+us+"B1"+us+mod_str;
	fname4 = fname4+us+argstr+us+"om"+us+om_str+us+"B1"+us+mod_str;
	printf("bimetric\n");
		

	fname = fname+extstr;
	fname2 = fname2+extstr;
	fname3 = fname3+extstr;
	fname4 = fname4+extstr;

	printf("%s\n",fname.c_str());
	printf("%s\n",fname2.c_str());
	printf("%s\n",fname3.c_str());
	printf("%s\n",fname4.c_str());

	FILE *fppass = fopen(fname.c_str(),"w");
	FILE *fppass2 = fopen(fname2.c_str(),"w");
	FILE *fppass3 = fopen(fname3.c_str(),"w");
	FILE *fp = fopen(fname4.c_str(),"w");

	cosmo_model_bigravity.run_cosmo(fppass,fppass2,fppass3);

	fclose(fppass);
	fclose(fppass2);
	fclose(fppass3);

	printf("RUN Done\n\n");

	for(theta = thetai;theta<=thetaend;theta*=dtheta)
	{
		
		theta_r = (M_PI/180.0)*theta;	
		//wgt = int_k(argc,argv,theta_r,cosmo_model_bigravity);
		wgt = int_logk(argc,argv,theta_r,cosmo_model_bigravity);
		
		fprintf(fp,"%lf\t%.10lf\n",theta,multi_fac*wgt);
		printf("%lf\t%.10lf\n",theta,multi_fac*wgt*1000000.0);
		
		

	}
	

	



	fclose(fp);


 test_pk(cosmo_model_bigravity.omega_dm_0,cosmo_model_bigravity.D1_0);

}








///##########################################################	
