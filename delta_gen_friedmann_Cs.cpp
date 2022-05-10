using namespace std;

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#define c_box 2.99

class cosmo_lcdm
{

	private:
	double omega_dm_0, H0, ratio;
	int model;

	public:
	cosmo_lcdm(double omega_dm_0_val=0.3,double h_val=0.7,int model=0)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = (h_val/c_box)*0.001;
		model = model;
		ratio = (1.0-omega_dm_0)/(omega_dm_0);
		printf("H0 %lf\n",H0);

	}

	double Hsqr(double a,double a0=1.0)
	{
		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		val = H0*H0*omega*( 1.0 + ratio*pow(a/a0,3.0) );
		return(val);


	}


	double Hb_t_by_Hb2(double a,double a0=1.0)
	{

		double val;
		
		val = -1.5/( 1.0 + ratio*pow(a/a0,3.0) );
		return(val);



	}

	double H_Diff(double a, double delta, double a0=1.0)
	{
		double diff;
		diff = -0.5*delta/( 1.0 + ratio*pow(a/a0,3.0) ) ;
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

	void pert_delta_aa(double *acc, double D[3], double D_a[3],double a,double a0=1.0)
	{
		double HbtbyHb2;	
		double c1,c2;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		
		
		c1 = -0.5/( 1.0 + ratio*pow(a/a0,3.0) );
		c2 = 0.0;

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);
		acc[2] = delta_aa(a, D[2],  D_a[2]);


	}

};


class cosmo_bigravity
{

	private:
	double omega_dm_0, H0, ratio,b,beta,gamma,B0;
	int model;

	public:
	double  B1;
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



class cosmo_dgp
{

	private:
	double omega_dm_0, H0, ratio, omega_r,beta, b;
	int model;

	public:
	cosmo_dgp(double omega_dm_0_val=0.3,double h_val=0.7,int model=1)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = (h_val/c_box)*0.001;
		model = model;
		omega_r = 0.5*(1.0-omega_dm_0)*0.5*(1.0-omega_dm_0);
		printf("Omega_r is %lf\n",omega_r);
		

	}

	double Hsqr(double a,double a0=1.0)
	{
		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		beta = sqrt( omega + omega_r);
		val = H0*H0*( omega + 2.0*sqrt(omega_r)*(sqrt(omega_r)+beta) );
		return(val);


	}


	double Hb_t_by_Hb2(double a,double a0=1.0)
	{

		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		beta = sqrt( omega + omega_r);
		
		val = -1.5*omega*(1.0+sqrt(omega_r)/beta)/( omega + 2.0*sqrt(omega_r)*(sqrt(omega_r)+beta) );
		return(val);



	}

	double H_Diff(double a, double delta, double a0=1.0)
	{
		double diff,omega,term1,term2;
		omega = omega_dm_0*pow(a0/a,3.0);
		beta = sqrt( omega + omega_r);
		b = sqrt(omega_r + omega*(1.0+delta));
		term1 = omega*(1.5*( sqrt(omega_r)*(1.0/beta -(1.0+delta)/b) ) - 0.5*delta );
		term2 = 2.0*sqrt(omega_r)*(b-beta);
		diff = (term1+term2)/( omega + 2.0*sqrt(omega_r)*(sqrt(omega_r)+beta) );
		return(diff);
		

	}

	double delta_aa(double a, double delta, double delta_a)
	{
		double HbtbyHb2, diff, acc,Hsqr_val,ddelta_dtau_sqr_by_adotsqr;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		diff = H_Diff(a, delta);
		Hsqr_val = Hsqr(a);
		//ddelta_dtau_sqr_by_adotsqr = delta_a*delta_a;

		acc = -(1.0/a)*(3.0 + HbtbyHb2)*delta_a + (4.0/3.0)*delta_a*delta_a/(1.0+delta) - 3.0*(1.0+delta)*diff/(a*a);
	
		return(acc);

	}


	double lin_delta_aa(double a, double delta, double delta_a,double a0=1.0)
	{
		double HbtbyHb2, diff, acc,theta,numer,denom;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		theta = pow(a0/a,3.0);
		numer = -theta*omega_dm_0*( 2.0*pow(omega_r,1.5) - sqrt(omega_r)*theta*omega_dm_0 + 2.0*omega_r*sqrt(omega_r + theta*omega_dm_0) 
						+ 2.0*theta*omega_dm_0*sqrt(omega_r + theta*omega_dm_0) );
		denom = 4.0*pow(omega_r+theta*omega_dm_0,1.5)*(2.0*omega_r + theta*omega_dm_0 + 2.0*sqrt(omega_r)*sqrt(omega_r+theta*omega_dm_0));
		diff = numer*delta/denom;
		acc = -(1.0/a)*(3.0 + HbtbyHb2)*delta_a - 3.0*diff/(a*a);
	
		return(acc);

	}


	void pert_delta_aa(double *acc, double D[3], double D_a[3],double a,double a0=1.0)
	{
		double HbtbyHb2, diff, theta,kappa1,lambda1,kappa2,lambda2,omega_dm;	
		double c1,c2;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		
		
		theta = pow(a0/a,3.0);
		omega_dm = omega_dm_0*theta;
		
		kappa1 = -( 2.0*(omega_dm_0 + omega_r/theta)*(omega_dm_0 + 4.0*omega_r/theta) +
					sqrt(omega_r)*sqrt( omega_r + omega_dm )*(5.0*omega_dm_0 + 8.0*omega_r/theta)/theta );
		
		lambda1 = 4.0*(omega_dm_0 + omega_r/theta)*(omega_dm_0 + omega_r/theta);
		
		
		c1 = kappa1/lambda1;

		kappa2 = omega_dm_0*omega_dm_0*sqrt(omega_r)*sqrt(omega_r + omega_dm)*(8.0*omega_r/theta - omega_dm_0);

	
		lambda2 = 16.0*(sqrt(omega_r) + sqrt( omega_r + omega_dm ))*(sqrt(omega_r) + sqrt( omega_r + omega_dm ))*
						(omega_dm_0 + omega_r/theta)*(omega_dm_0 + omega_r/theta)*(omega_dm_0 + omega_r/theta);

		c2 = kappa2/lambda2;

		diff = H_Diff(a, D[2]);

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);
		acc[2] = delta_aa(a, D[2],  D_a[2]);
			 
			  



	}


};



int f_sigma_cal(int argc,char *argv[],double sig8,double *x_lst,double *fs8,int n,int ax=1)
{
	double D[3],D_a[3],D_rk[5][3], D_a_rk[5][3], acc[3];
	double D_i[3], D_a_i[3];

	double a,ai,ai_burn,a0,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};
	double om_dm_0, model_param;
	double D1_0;
	double a_lst[n];	
	int a_lst_cntr=0;
	double a_lst_now ;
	

	int i,j,theory;

	for(i=0;i<n;++i)
	{
		if(ax)
		a_lst[i] = x_lst[i];
		
		else
		a_lst[i] = 1.0/(1.0+x_lst[i]);


	}


	string fname = "delta";
	string argstr = argv[1];
	string om_str = argv[2];
	string mod_str = argv[3];
	string extstr = ".txt";
	string us = "_";
	
	om_dm_0 = atof(argv[2]);
	model_param = atof(argv[3]);
	

	
	cosmo_lcdm cosmo_model_lcdm(om_dm_0);
	cosmo_dgp cosmo_model_dgp(om_dm_0);
	cosmo_bigravity cosmo_model_bigravity(om_dm_0,model_param);
	
	if(!strcmp(argv[1],"lcdm"))
	{	fname = fname+us+argstr+us+"om"+us+om_str;
		printf("lcdm\n");
		theory = 0;
	}
	else
	if(!strcmp(argv[1],"dgp"))
	{	
		fname = fname+us+argstr+us+"om"+us+om_str;
		printf("dgp\n");
		theory = 1;
	}
	else
	if(!strcmp(argv[1],"bimetric"))
	{	
		fname = fname+us+argstr+us+"om"+us+om_str+us+"B1"+us+mod_str;
		printf("bimetric\n");
		
		theory = 2;
	}

	string fname_fs8 = fname+"_fs8";
	fname = fname+extstr;
	fname_fs8 = fname_fs8+extstr;

	FILE *fp = fopen(fname.c_str(),"w");
	FILE *fp_fs8 = fopen(fname_fs8.c_str(),"w");
	printf("%s %s\n",fname.c_str(),fname_fs8.c_str());

	da = 0.0000001;
	ai = 0.001;
	ai_burn = ai*0.01;
	a0 = 1.0;


	D_a_i[0] = 1.0;
	D_i[0] = ai_burn*D_a_i[0];
	
	D[0] = D_i[0];
	D_a[0] = D_a_i[0];

	D_i[1] = (34.0/7.0)*D_i[0]*D_i[0]/3.0;
	D_a_i[1] = D_i[1]/ai;
	
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

		if(a>=a_lst_now)
		{
			fs8[a_lst_cntr] = sig8*a*D_a[0];

			++a_lst_cntr;
			a_lst_now = a_lst[a_lst_cntr];
		}


		if((a>ai)&&(burn))
		{
			burn = 0;

			D_i[0] = D[0];
			D_a_i[0] = D_a[0];
	
			D_i[2] = D[2];
			D_a_i[2] = D_a[2];

			//D[1] = 4.858*D[0]*D[0]/3.0;
			//D_a[1] = D_a[0];//1.0;
		
			D_i[1] = D[1];
			D_a_i[1] = D_a[1];
		
			printf("intial dc is %lf  %lf  %lf\n",D_i[0],D_i[1],D_i[2]);


		}



		D_rk[0][0] = D[0];
		D_a_rk[0][0] = D_a[0]; 

		D_rk[0][1] = D[1];
		D_a_rk[0][1] = D_a[1]; 

		D_rk[0][2] = D[2];
		D_a_rk[0][2] = D_a[2]; 
		
		if(!(burn)&&((cntr%100)==0))
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			a,D[0],D[1],D[2],D[0]*ai/(a*D_i[0]),D[1]*ai/(a*D_i[1]),D[0]/D_i[0],D[1]/D_i[1],D[2]/D_i[2],
			3.0*D[1]/(D[0]*D[0]),3.0*D[2]/(D[0]*D[0]));

		

		for(i=1;i<=4;++i)
		{
			if(theory==0)	
			cosmo_model_lcdm.pert_delta_aa(acc, D_rk[0], D_a_rk[0],a);		
			
			if(theory==1)			
			cosmo_model_dgp.pert_delta_aa(acc, D_rk[0], D_a_rk[0],a);

			if(theory==2)			
			cosmo_model_bigravity.pert_delta_aa(acc, D_rk[0], D_a_rk[0],a);

			
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

	D1_0 = D[0];

	for(i=0;i<n;++i)
	{

		fs8[i] = fs8[i]/D1_0;
		if(ax)
		fprintf(fp_fs8,"%lf\t%lf\n",a_lst[i],fs8[i]);
		else
		fprintf(fp_fs8,"%lf\t%lf\n",x_lst[i],fs8[i]);
	}


	fclose(fp);
	fclose(fp_fs8);
	
	return(1);


}


double * read_fs8_data(FILE *fp_data,int &n,int up_bnd = 100)
{

	int fpr=1,data_cnt=0;
	double *fs = new double[up_bnd];
	double *erl = new double[up_bnd];
	double *eru = new double[up_bnd];


	while(fpr>0)
	{
		fpr = fscanf(fp_data,"%lf\t%lf\t%lf\n",&fs[data_cnt],&erl[data_cnt],&eru[data_cnt]);
		if(fpr>0)
		++data_cnt;
		//if(data_cnt<10)
		//printf("%d %d \n",fpr,data_cnt);
		

	}


	double *data = new double[data_cnt];
	for(int i=0;i<data_cnt;++i)
	{
		data[i] = fs[i];
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
	double z_l[12] = {0.02,0.08,0.15,0.19,0.23,0.35,0.45,0.55,0.65,0.78,0.89,0.999};
	int n = 12,b;
	int ax = 0;
	double *data;

	double fs8[12];

	FILE *fpdata = fopen("data.txt","r");
	data = read_fs8_data(fpdata,n);
	printf("n is %d\n",n);
	for(int i=0;i<n;++i)
	{
		
		printf("%d %lf\n",i,data[i]);
		

	}

	

	b = f_sigma_cal( argc,argv, sig8,a_l,fs8,n,ax);

}








///##########################################################	
