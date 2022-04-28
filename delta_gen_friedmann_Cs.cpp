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

	void pert_delta_aa(double *acc, double D[2], double D_a[2],double a,double a0=1.0)
	{
		double HbtbyHb2;	
		double c1,c2;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		
		
		c1 = -0.5/( 1.0 + ratio*pow(a/a0,3.0) );
		c2 = 0.0;

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);



	}

};


class cosmo_bigravity
{

	private:
	double omega_dm_0, H0, ratio,b,beta,gamma,B0,B1;
	int model;

	public:
	cosmo_bigravity(double omega_dm_0_val=0.3,double h_val=0.7,int model=2)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = (h_val/c_box)*0.001;
		model = model;
		B1 = 0.5;
		B0 = 3.0*(1.0-omega_dm_0-B1*B1/3.0);
		ratio = (1.0-omega_dm_0)/(omega_dm_0);
		printf("H0 %lf\n",H0);

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


	void pert_delta_aa(double *acc, double D[2], double D_a[2],double a,double a0=1.0)
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
		lambda2 = lambda2*(3.0 - B1*B1 - 3.0*theta + 3.0*omega_dm_0*theta + 6.0*sqrt(B1*B1/3.0 + (B0/6.0 + 0.5*omega_dm_0*theta)*(B0/6.0 + 0.5*omega_dm_0*theta) ));
	
		c2 = kappa2/lambda2;

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);



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


	void pert_delta_aa(double *acc, double D[2], double D_a[2],double a,double a0=1.0)
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

		acc[0] = -(3.0 + HbtbyHb2)*D_a[0]/a - 3.0*c1*D[0]/(a*a);
		acc[1] = -(3.0 + HbtbyHb2)*D_a[1]/a + (8.0/3.0)*D_a[0]*D_a[0] - 3.0*c1*D[1]/(a*a) - 6.0*(c1+c2)*D[0]*D[0]/(a*a);



	}


};







int main(int argc,char *argv[])
{

	double D[2],D_a[2],D_rk[2][5], D_a_rk[2][5], acc[2];
	double D_i[2], D_a_i[2];

	double a,ai,ai_burn,a0,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};

	int i,j,theory;

	string fname = "delta_";
	string argstr = argv[1];
	string extstr = "_B1_0p5.txt";
	fname = fname+argstr+extstr;
	printf("%s\n",fname.c_str());

	FILE *fp = fopen(fname.c_str(),"w");
	cosmo_lcdm cosmo_model_lcdm(1.0);
	cosmo_dgp cosmo_model_dgp(0.3);
	cosmo_bigravity cosmo_model_bigravity(0.3);
	
	if(!strcmp(argv[1],"lcdm"))
	{
		printf("lcdm\n");
		theory = 0;
	}
	else
	if(!strcmp(argv[1],"dgp"))
	{	printf("dgp\n");
		theory = 1;
	}
	else
	if(!strcmp(argv[1],"bimetric"))
	{	printf("bimetric\n");
		theory = 2;
	}

	da = 0.0000001;
	ai = 0.001;
	ai_burn = ai*0.01;
	a0 = 1.0;


	D_a_i[0] = 0.01;
	D_i[0] = ai_burn*D_a_i[0];
	
	D[0] = D_i[0];
	D_a[0] = D_a_i[0];

	D_i[1] = (34.0/7.0)*D_i[0]*D_i[0]/3.0;
	D_a_i[1] = D_i[1]/ai;
	
	D[1] = D_i[1];
	D_a[1] = D_a_i[1];



	int burn = 1, cntr;

	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);

	for(a=ai_burn,cntr=0;a<=a0;a+=da,++cntr)
	{
		ak = a;


		if((a>ai)&&(burn))
		{
			burn = 0;

			D_i[0] = D[0];
			D_a_i[0] = D_a[0];

			//D[1] = 4.858*D[0]*D[0]/3.0;
			//D_a[1] = D_a[0];//1.0;
		
			D_i[1] = D[1];
			D_a_i[1] = D_a[1];
		
			printf("intial dc is %lf  %lf\n",D_i[0],D_i[1]);


		}



		D_rk[0][0] = D[0];
		D_a_rk[0][0] = D_a[0]; 

		D_rk[1][0] = D[1];
		D_a_rk[1][0] = D_a[1]; 
		
		if(!(burn)&&((cntr%100)==0))
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,D[0],D[1],D[0]*ai/(a*D_i[0]),D[1]*ai/(a*D_i[1]),D[0]/D_i[0],D[1]/D_i[1],3.0*D[1]/(D[0]*D[0]));

		

		for(i=1;i<=4;++i)
		{
			if(theory==0)	
			cosmo_model_lcdm.pert_delta_aa(acc, D, D_a,a);		
			
			if(theory==1)			
			cosmo_model_dgp.pert_delta_aa(acc, D, D_a,a);

			if(theory==2)			
			cosmo_model_bigravity.pert_delta_aa(acc, D, D_a,a);

			
		   for(j=0;j<2;++j)
			{D_a_rk[j][i] = da*acc[j];
			 D_rk[j][i] = da*D_a_rk[j][0];
				

			 

			 D_rk[j][0] = D[j] + rk_coef[i-1]*D_rk[j][i];
			 D_a_rk[j][0] = D_a[j] + rk_coef[i-1]*D_a_rk[j][i];
			}

			ak = a + rk_coef[i-1]*da;
			if(a==ai)
			 printf("acc %.10lf  %.10lf\n",D_rk[0][i],D_rk[1][i]);

		}


	    for(j=0;j<2;++j)
		{
		 D[j] = D[j] + (1.0/6.0)*(D_rk[j][1]+2.0*D_rk[j][2]+2.0*D_rk[j][3]+D_rk[j][4]);
		 D_a[j] = D_a[j] + (1.0/6.0)*(D_a_rk[j][1]+2.0*D_a_rk[j][2]+2.0*D_a_rk[j][3]+D_a_rk[j][4]);

		}


	}


	

}
