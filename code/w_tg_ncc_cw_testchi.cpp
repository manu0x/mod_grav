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










class cosmo_negcc_cw
{

	private:
	double  ratio_dem,ratio_deL,b,beta,w;
	int model;

	public:
	double omega_dm_0, H0,omega_L_0,omega_de_0;
	double D1_0,net_chi,chi0ret; 
	
	double dai;
	

	double *cs_D1,*cs_f,*cs_cmd,*a_sp,*D1_sp,*f_sp,*cmd_sp;
	cosmo_negcc_cw(double omega_dm_0_val=0.3,double olparam=0.0,double wparam=-1.0,double h_val=0.7,int model=2)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = (h_val/c_box)*0.001;
		model = model;
		w = wparam;
		omega_L_0 = olparam;
		omega_de_0 = (1.0-omega_dm_0-omega_L_0);
		ratio_dem = omega_de_0/(omega_dm_0);
		ratio_deL = omega_L_0/(omega_dm_0);

		chi0ret = (2.0/3.0)*atanh(sqrt(1.0-omega_dm_0))/sqrt(1.0-omega_dm_0);
		printf("omL %lf  w  %lf  ratio_dem %lf \n",omega_L_0,w,ratio_dem );

		

		TFmdm_set_cosm(omega_dm_0, 0.0223/(h_val*h_val), 0.0,
		1, (1.0-omega_dm_0), h_val, 1.0);

	}

	double Hsqr(double a,double a0=1.0)
	{
		double val,omega_dm,omega_de;
		omega_dm = omega_dm_0*pow(a0/a,3.0);
		omega_de = omega_de_0*pow(a0/a,3.0*(1.0+w));
		val = H0*H0*( omega_dm + omega_de + omega_L_0);
		return(val);


	}


	void g(double x,double dr[4])
	{
 	   
 	   
    	   dr[0] = x + omega_de_0*pow(x/omega_dm_0,1.0+w) + omega_L_0  ;
    	   dr[1] = 1.0 + (1.0+w)*ratio_dem*pow(x/omega_dm_0,w);
    	   dr[2] = (1.0+w)*w*ratio_dem*pow(x/omega_dm_0,w-1.0)/omega_dm_0;
    	   dr[3] = (1.0+w)*w*(w-1.0)*ratio_dem*pow(x/omega_dm_0,w-2.0)/(omega_dm_0*omega_dm_0);



	}

	double lcdm_chi(double a,double a0=1.0)
	{

		double targ,ret;
		targ = sqrt((1.0-omega_dm_0)/((1.0-omega_dm_0) + omega_dm_0*pow(a0/a,3.0)));

		ret = chi0ret - (2.0/3.0)*atanh(targ)/sqrt(1.0-omega_dm_0);
		return(ret/H0byc);
	}



	double Hb_t_by_Hb2(double a,double a0=1.0)
	{

		double val;
		
		val = -1.5*(1.0+ratio_dem*pow(a,-3.0*w)*(1.0+w))/( 1.0 + ratio_dem*pow(a,-3.0*w) + ratio_deL*pow(a,3.0) );
		return(val);



	}

	double H_Diff(double a, double delta, double a0=1.0)
	{
		double diff;
		diff = -0.5*delta/( 1.0 + ratio_dem*pow(a,-3.0*w) + ratio_deL*pow(a,3.0)  ) ;
		//printf("lcdm diff %lf\n",diff);
		return(diff);
		

	}

	double delta_aa(double a, double delta, double delta_a)
	{
		double HbtbyHb2, diff, acc,Hsqr_val,ddelta_dtau_sqr_by_adotsqr;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		diff = H_Diff(a, delta);
		Hsqr_val = Hsqr(a);
		//ddelta_dtau_sqr_by_adotsqr = a*a*delta_a*delta_a;

		acc = -(1.0/a)*(3.0 + HbtbyHb2)*delta_a + (4.0/3.0)*delta_a*delta_a/(1.0+delta) - 3.0*(1.0+delta)*diff/(a*a);
	
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
		double HbtbyHb2;	
		double c1,c2;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		
		
		c1 = -0.5/( 1.0 + ratio_dem*pow(a,-3.0*w) + ratio_deL*pow(a,3.0)  );
		c2 = 0.0;

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a); //printf("acc[0] %.10lf\t%.10lf\t%.10lf\n",acc[0],D_a[0]/a,D[0]/(a*a));
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);
		acc[2] = delta_aa(a, D[2],  D_a[2]);
		acc[3] = delta_lin_lue(a,D[3],D_a[3],1);
		acc[4] = delta_lin_lue(a,D[4],D_a[4],0);


	}


     int run_cosmo(FILE *fp,FILE *fp_de,FILE *fp_mg,int de_mode=1,double da=0.000001)
     {	printf("run cosmo omL %lf  w  %lf\n",omega_L_0,w);

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
	
		ak = a;
		
		//
		//if(((cntr%1000)==0))
		if(!(burn)&&((cntr%1000)==0))
		{	x =  omega_dm_0*pow(a0/a,3.0);
			g(x,gv);


			fv = a*(D_a[0]/D[0]);
			isw_potn = (1.0-fv)*D[0];	

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
	isw_potn = (1.0-fv)*D[0];	

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

double int_z(int argc,char *argv[],double theta,double k,cosmo_negcc_cw bimet,double da =  0.00001)
{


  int i,j,cur_i;
  double in_z,D1,chi,gv[4],fv,j0,wgv,ktc,x;
  double ai,a0,a,delta;
  double cur_int;
  double as[3];
  double chical;
 
  double dsmp[3]={0.0,0.5*da,da};
  double smp[3]={1.0,4.0,1.0};

  int cntr;

  ai = 0.001;
  a0 = 1.0;

  in_z = 0.0;

  FILE *fpchi = fopen("chitest.txt","w");
  
  for(a=ai,cntr=0;a<=a0;a+=da,++cntr)
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

		chical = bimet.lcdm_chi(as[j]);
		if((j==0)&&((cntr%1000)==0))
		fprintf(fpchi,"%lf\t%lf\t%lf\n",as[j],bimet.net_chi-chi,chical);
		j0 = gsl_sf_bessel_J0(ktc);
		wgv = wg(a0/as[j]-1.0);

		cur_int = sqrt(gv[0])*(D1/as[j])*(D1/as[j])*( (1.0-fv))*wgv*j0;

		in_z+=(da/6.0)*smp[j]*cur_int;
			
		

	}

	


  }

	fclose(fpchi);
  

  return(in_z);



}






double int_k(int argc,char *argv[],double theta,cosmo_negcc_cw bimet,double kini=0.00004,double kend=3.0,double dk=0.01)
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

	   in_k1 = int_z(argc,argv,theta,k,bimet);//printf("theta %lf  k_frac %lf  %lf\n",theta,k/kend,in_k1);
	   in_k1 = pk1*in_k1/k;

	   
	
	   in_k2 = int_z(argc,argv,theta,k+0.5*dk,bimet);
	   in_k2 = pk2*in_k2/(k+0.5*dk);

	   in_k3 = int_z(argc,argv,theta,k+dk,bimet);
	   in_k3 = pk3*in_k3/(k+dk);

	   in_k+= (dk/6.0)*(in_k1+4.0*in_k2+in_k3);


	}
	
	return(in_k);


}


double int_logk(int argc,char *argv[],double theta,cosmo_negcc_cw bimet,double lkini=-12.0,double lkend=1.0,double dlk=0.05)
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


	double om_dm_0, model_param,model_param2;
	double D1_0; int theory;

	double multi_fac,T0,bias;
	double theta,thetai,thetaend,dtheta,theta_r,wgt;

	// These thetas are in degrees while for calc. we convert value to rad in loop before passing rad val to intgrat fncts
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
	string mod_str2 = argv[4];
	string extstr = ".txt";
	string us = "_";

	
	
	om_dm_0 = atof(argv[2]);
	model_param = atof(argv[3]);
	model_param2 = atof(argv[4]);

	T0 = 2.72548;
	bias = 5.47;

	multi_fac = 3.0*T0*twopie*twopie*bias*om_dm_0*H0byc*H0byc*H0byc;
	//multi_fac = 3.0*T0*bias*om_dm_0*H0byc*H0byc*H0byc/twopie;
	

	

	cosmo_negcc_cw cosmo_model_nccw(om_dm_0,model_param,model_param2);

	
	

	

	fname = fname+us+argstr+us+"om"+us+om_str+us+"omL"+us+mod_str+us+"w"+us+mod_str2;
	fname2 = fname2+us+argstr+us+"om"+us+om_str+us+"omL"+us+mod_str+us+"w"+us+mod_str2;
	fname3 = fname3+us+argstr+us+"om"+us+om_str+us+"omL"+us+mod_str+us+"w"+us+mod_str2;
	fname4 = fname4+us+argstr+us+"om"+us+om_str+us+"omL"+us+mod_str+us+"w"+us+mod_str2;

	printf("-CC constant w\n");
		

	fname = fname+extstr;
	fname2 = fname2+extstr;
	fname3 = fname3+extstr;
	fname4 = fname4+extstr;

	printf("%s\n",fname.c_str());
	printf("%s\n",fname2.c_str());

	printf("\nomdm0 %lf\n",om_dm_0);
	printf("omL0 %lf\n",model_param);
	printf("w %lf\n\n",model_param2);

	FILE *fppass = fopen(fname.c_str(),"w");
	FILE *fppass2 = fopen(fname2.c_str(),"w");
	FILE *fppass3 = fopen(fname3.c_str(),"w");
	FILE *fp = fopen(fname4.c_str(),"w");

	cosmo_model_nccw.run_cosmo(fppass,fppass2,fppass3);

	fclose(fppass);
	fclose(fppass2);
	fclose(fppass3);

	printf("RUN Done\n\n");

/*	for(theta = thetai;theta<=thetaend;theta*=dtheta)
	{
		theta_r = (M_PI/180.0)*theta;				//Converting theta into radians
		//wgt = int_k(argc,argv,theta_r,cosmo_model_nccw);
		wgt = int_logk(argc,argv,theta_r,cosmo_model_nccw);
		
		fprintf(fp,"%lf\t%.10lf\n",theta,multi_fac*wgt);
		printf("%.10lf\t%.10lf\n",theta,multi_fac*wgt*1000000.0);
		

	}
	
*/
	
	wgt = int_logk(argc,argv,(M_PI/180.0)*0.1,cosmo_model_nccw);
		
		
		printf("%.10lf\t%.10lf\n",0.1,multi_fac*wgt*1000000.0);

	
	fclose(fp);


 test_pk(cosmo_model_nccw.omega_dm_0,cosmo_model_nccw.D1_0);

}













