using namespace std;

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <gsl/gsl_sf_bessel.h>

#define c_box 2.99
#define twopie 2.0*M_PI




class cosmo_bigravity
{

	private:
	double  ratio,b,beta,gamma;
	int model;

	public:
	double omega_dm_0, H0,B1;
	double  B0;
	cosmo_bigravity(double omega_dm_0_val=0.3,double B1param=0.0,double h_val=0.7,int model=2)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = (h_val/c_box)*0.001;
		model = model;
		B1 = B1param;
		B0 = 3.0*(1.0-omega_dm_0-B1*B1/3.0);
		ratio = (1.0-omega_dm_0)/(omega_dm_0);
		printf("H0 %lf  B1  %lf\n",H0,B1);

	}

	double Hsqr(double a,double a0=1.0)
	{
		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		val = H0*H0*(   0.5*omega + B0/6.0  + sqrt((0.5*omega + B0/6.0)*(0.5*omega + B0/6.0) + B1*B1/3.0 ) );
		return(val);


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


	void pert_delta_aa(double *acc, double D[3], double D_a[3],double a,double a0=1.0)
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
		lambda2 = lambda2*(3.0 - B1*B1 - 3.0*omega_dm_0 + 3.0*omega_dm_0*theta + 6.0*sqrt(B1*B1/3.0 + (B0/6.0 + 0.5*omega_dm_0*theta)*(B0/6.0 + 0.5*omega_dm_0*theta) ));
	
		c2 = kappa2/lambda2;

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);
		acc[2] = delta_aa(a, D[2],  D_a[2]);



	}



	
	

};





void g(double x,double dr[4],double B0,double B1=0.0)
{
    double b0t,b1t;
    b0t = B0/6.0 + x;
    b1t = B1*B1/3.0;
    dr[0] =  b0t + sqrt(b0t*b0t + b1t);
    dr[1] = 0.5 + b0t/sqrt(b1t + b0t*b0t);
    dr[2] = -b0t*b0t/pow(b1t + b0t*b0t,1.5) + 1.0/sqrt(b1t + b0t*b0t);
    dr[3] = 3.0*b0t*b0t*b0t/pow(b1t + b0t*b0t,2.5) -3.0*b0t/pow(b1t + b0t*b0t,1.5);	



}


double wg(double z)
{
	return(1.0);
}

double Pk(double k)
{

	return (1.0/k);

}





double int_z(int argc,char *argv[],double theta,double k,cosmo_bigravity bimet,FILE *fp)
{
	double D[3],D_a[3],D_rk[5][3], D_a_rk[5][3], acc[3];
	double D_i[3], D_a_i[3];

	double a,ai,ai_burn,a0,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};
	

	
	double in_z,chi1,chi2,chi3,chi,cur_z,gv[4],x,wgv,ktc,j0,fv;
	

	int i,j,theory,acntr;

	



	da = 0.000005;
	ai = 0.001;
	ai_burn = ai*0.1;
	a0 = 1.0;


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
	
	D[2] = D_i[2];
	D_a[2] = D_a_i[2];



	int burn = 1, cntr;


	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);

	for(a=ai_burn,cntr=0;a<=a0;a+=da,++cntr)
	{
		ak = a;




		if((a>ai)&&(burn))
		{
			burn = 0;
			acntr = 0;

			D_i[0] = D[0];
			D_a_i[0] = D_a[0];
	
			D_i[2] = D[2];
			D_a_i[2] = D_a[2];

			//D[1] = 4.858*D[0]*D[0]/3.0;
			//D_a[1] = 4.858*2.0*D_a[0]*D[0]/3.0;
		
			//D_i[1] = D[1];
			//D_a_i[1] = D_a[1];
		
			//printf("intial dc is %lf  %lf  %lf\n",D_i[0],D_i[1],D_i[2]);
			in_z = 0.0;
			chi = 0.0;
			


		}



		D_rk[0][0] = D[0];
		D_a_rk[0][0] = D_a[0]; 

		D_rk[0][1] = D[1];
		D_a_rk[0][1] = D_a[1]; 

		D_rk[0][2] = D[2];
		D_a_rk[0][2] = D_a[2]; 
		
		//
		//if(((cntr%1000)==0))
		if(!(burn)&&((cntr%1000)==0))
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			a,D[0],D[1],D[2],D[0]*ai/(a*D_i[0]),D[1]*ai/(a*D_i[1]),D[0]/D_i[0],D[1]/D_i[1],D[2]/D_i[2],
			3.0*D[1]/(D[0]*D[0]),3.0*D[2]/(D[0]*D[0]));

		if(a>=ai)
		{
		  x =  bimet.omega_dm_0*pow(a0/a,3.0);
		  chi1 = 1.0/(a*a*sqrt(bimet.Hsqr(a)));
		  chi2 = 1.0/((a+0.5*da)*(a+0.5*da)*sqrt(bimet.Hsqr(a+0.5*da)));
		  chi3 = 1.0/((a+da)*(a+da)*sqrt(bimet.Hsqr(a+da)));
		  chi+= (da/6.0)*(chi1+4.0*chi2+chi3)*(bimet.H0);
		  g(x,gv, bimet.B0, bimet.B1);
		  fv = a*(D_a[0]/D[0]);
		  ktc = k*theta*chi;
		  j0 = gsl_sf_bessel_J0(ktc);
		  wgv = wg(a0/a-1.0);
		  
		  cur_z = sqrt(gv[0])*D[0]*D[0]*( (1.0-fv)*(gv[1] + 1.5*x*gv[2]) + 1.5*x*(5.0*gv[2] + 3.0*x*gv[3]) )*wgv*j0;

		  if(acntr%2==0)
		   in_z+=da*2.0*cur_z/3.0;
		  else
		   in_z+=da*4.0*cur_z/3.0;

		  if((acntr == 0)||(a>=a0))
			in_z+= da*cur_z/3.0;

		}
		

		for(i=1;i<=4;++i)
		{
					
			bimet.pert_delta_aa(acc, D_rk[0], D_a_rk[0],a);

			
		   for(j=0;j<3;++j)
			{D_a_rk[i][j] = da*acc[j];
			 D_rk[i][j] = da*D_a_rk[0][j];
				

			 

			 D_rk[0][j] = D[j] + rk_coef[i-1]*D_rk[i][j];
			 D_a_rk[0][j] = D_a[j] + rk_coef[i-1]*D_a_rk[i][j];
			}

			ak = a + rk_coef[i-1]*da;
			if(a==ai)
			 printf("acc %.10lf  %.10lf\n",D_rk[i][0],D_rk[i][1]);

		}


	    for(j=0;j<3;++j)
		{
		 D[j] = D[j] + (1.0/6.0)*(D_rk[1][j]+2.0*D_rk[2][j]+2.0*D_rk[3][j]+D_rk[4][j]);
		 D_a[j] = D_a[j] + (1.0/6.0)*(D_a_rk[1][j]+2.0*D_a_rk[2][j]+2.0*D_a_rk[3][j]+D_a_rk[4][j]);

		

		}


	}

  return(in_z);



}







double int_k(int argc,char *argv[],double theta,cosmo_bigravity bimet,FILE *fp,double kini=1.0,double kend=2.0,double dk=0.01)
{

	double in_k,in_k1,in_k2,in_k3;
	double pk1,pk2,pk3;
	double k;
	in_k = 0.0;
	
	printf("\n\nTHETA %lf\n",theta);
	
	for(k=kini;k<=kend;k+=dk)
	{ 

	   pk1 = Pk(k);
	   pk2 = Pk(k+0.5*dk);
	   pk3 = Pk(k+dk);

	   in_k1 = int_z(argc,argv,theta,k,bimet,fp);printf("theta %lf  k_frac %lf  %lf\n",theta,k,in_k1);
	   in_k1 = pk1*in_k1/k;
	
	   in_k2 = int_z(argc,argv,theta,k+0.5*dk,bimet,fp);
	   in_k2 = pk2*in_k2/(k+0.5*dk);

	   in_k3 = int_z(argc,argv,theta,k+dk,bimet,fp);
	   in_k3 = pk3*in_k3/(k+dk);

	   in_k+= (dk/6.0)*(in_k1+4.0*in_k2+in_k3);


	}
	
	return(in_k);


}





//int f_sigma_cal(int argc,char *argv[],double sig8,double *al,double *fs8,int n,int ax=1)



double * read_fs8_data(FILE *fp_data,int &n,int up_bnd = 100)
{

	int fpr=1,data_cnt=0;
	double *fs = new double[up_bnd];
	double *xs = new double[up_bnd];
	double *erl = new double[up_bnd];
	double *eru = new double[up_bnd];


	while(fpr>0)
	{
		fpr = fscanf(fp_data,"%lf\t%lf\t%lf\n",&xs[data_cnt],&fs[data_cnt],&erl[data_cnt],&eru[data_cnt]);
		if(fpr>0)
		++data_cnt;
		//if(data_cnt<10)
		//printf("%d %d \n",fpr,data_cnt);
		

	}


	double *data = new double[data_cnt*2];
	
	for(int i=0,j=0;i<data_cnt;++i,j+=2)
	{
		data[j] = xs[i];
		data[j+1] = fs[i];
		//printf("%d %lf %lf %lf\n",i,data[i],erl[i],eru[i]);
		

	}
	
	n = data_cnt;

	delete[] fs;
	delete[] erl;
	delete[] eru;

	return data;

}


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
	double theta,thetai,thetaend,dtheta,wgt;

	
	thetai = 0.1;
	thetaend = 0.1;
	dtheta = 0.25;

	string fname = "wgtxt";
	string argstr = argv[1];
	string om_str = argv[2];
	string mod_str = argv[3];
	string extstr = ".txt";
	string us = "_";

	
	
	om_dm_0 = atof(argv[2]);
	model_param = atof(argv[3]);

	T0 = 2.14;
	bias = 5.47;

	multi_fac = 3.0*T0*twopie*twopie*bias*om_dm_0;
	

	

	cosmo_bigravity cosmo_model_bigravity(om_dm_0,model_param);

	
	

	

	fname = fname+us+argstr+us+"om"+us+om_str+us+"B1"+us+mod_str;
	printf("bimetric\n");
		

	fname = fname+extstr;

	printf("%s\n",fname.c_str());

	FILE *fppass = fopen(fname.c_str(),"w");
	FILE *fp = fopen("wgt_lcdm2.txt","w");

	for(theta = thetai;theta<=thetaend;theta+=dtheta)
	{
		
		
		wgt = int_k(argc,argv,theta,cosmo_model_bigravity,fppass);
		
		fprintf(fp,"%lf\t%lf\n",theta,wgt);
		printf("%lf\t%lf\n",theta,wgt);
		

	}
	
	

	fclose(fppass);
	fclose(fp);
}








///##########################################################	
