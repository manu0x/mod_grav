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

	double lin_delta_aa(double a, double delta, double delta_a,double a0=1.0)
	{
		double HbtbyHb2, diff, acc,Hsqr_val,ddelta_dtau_sqr_by_adotsqr;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		diff = H_Diff(a, delta);
		Hsqr_val = Hsqr(a);
		//ddelta_dtau_sqr_by_adotsqr = a*a*delta_a*delta_a;

		acc = -(1.0/a)*(3.0 + HbtbyHb2)*delta_a - 3.0*diff/(a*a);
	
		return(acc);

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
		B1 = 0.00;
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


};







int main(int argc,char *argv[])
{

	double delta,delta_a,delta_rk[5], delta_a_rk[5], acc;
	double delta_i, delta_a_i;

	double a,ai,ai_burn,a0,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};

	int i,theory;

	string fname = "delta_";
	string argstr = argv[1];
	string extstr = ".txt";
	fname = fname+argstr+extstr;
	printf("%s\n",fname.c_str());

	FILE *fp = fopen(fname.c_str(),"w");
	cosmo_lcdm cosmo_model_lcdm(0.3);
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
	ai_burn = 0.1*ai;
	a0 = 1.0;

	delta_i = ai;
	delta_a_i = 0.000;
	
	delta = delta_i;
	delta_a = delta_a_i;



	int burn = 1, cntr;

	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);

	for(a=ai_burn,cntr=0;a<=a0;a+=da,++cntr)
	{
		ak = a;
		delta_rk[0] = delta;
		delta_a_rk[0] = delta_a; 
		
		if(!(burn)&&((cntr%100)==0))
		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta*ai/(a*delta_i),delta/delta_i,delta_a/delta_a_i);

		if((a>ai)&&(burn))
		{
			burn = 0;
			delta_i = delta;
			delta_a_i = delta_a;
			printf("intial dc is %lf\n",delta_i);


		}

		for(i=1;i<=4;++i)
		{
			if(theory==0)			
			acc  = cosmo_model_lcdm.lin_delta_aa( ak,  delta_rk[0], delta_a_rk[0]);
			if(theory==1)			
			acc  = cosmo_model_dgp.lin_delta_aa( ak,  delta_rk[0], delta_a_rk[0]);
			if(theory==2)			
			acc  = cosmo_model_bigravity.lin_delta_aa( ak,  delta_rk[0], delta_a_rk[0]);

			

			delta_a_rk[i] = da*acc;
			delta_rk[i] = da*delta_a_rk[0];

			if(a==ai)
			printf("acc %.10lf\n",delta_rk[i]);

			delta_rk[0] = delta + rk_coef[i-1]*delta_rk[i];
			delta_a_rk[0] = delta_a + rk_coef[i-1]*delta_a_rk[i];

			ak = a + rk_coef[i-1]*da;


		}

		delta = delta + (1.0/6.0)*(delta_rk[1]+2.0*delta_rk[2]+2.0*delta_rk[3]+delta_rk[4]);
		delta_a = delta_a + (1.0/6.0)*(delta_a_rk[1]+2.0*delta_a_rk[2]+2.0*delta_a_rk[3]+delta_a_rk[4]);

		


	}


	

}
