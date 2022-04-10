using namespace std;

#include <math.h>
#include <stdio.h>
#include <string.h>

class cosmo_lcdm
{

	private:
	double omega_dm_0, H0, ratio;
	int model;

	public:
	cosmo_lcdm(double omega_dm_0_val=0.3,double h_val=0.7,int model=0)
	{
		omega_dm_0 = omega_dm_0_val;
		//H0 = h_val;
		model = model;
		ratio = (1.0-omega_dm_0)/(omega_dm_0);

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


};


class cosmo_dgp
{

	private:
	double omega_dm_0, H0, ratio, omega_r,beta, b;
	int model;

	public:
	cosmo_dgp(double omega_dm_0_val=0.3,double h_val=0.7,int model=0)
	{
		omega_dm_0 = omega_dm_0_val;
		//H0 = h_val;
		model = model;
		omega_r = 0.5*(1.0-omega_dm_0)*0.5*(1.0-omega_dm_0);
		

	}


	double Hb_t_by_Hb2(double a,double a0=1.0)
	{

		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		beta = sqrt( omega + omega_r);
		
		val = -1.5*omega*(1.0+sqrt(omega_r)/beta)/( omega + 2.0*omega_r*(1.0+beta/sqrt(omega_r)) );
		return(val);



	}

	double H_Diff(double a, double delta, double a0=1.0)
	{
		double diff,omega,term1,term2;
		omega = omega_dm_0*pow(a0/a,3.0);
		beta = sqrt( omega + omega_r);
		b = sqrt(omega_r + omega*(1.0+delta));
		term1 = omega*(1.5*( sqrt(omega_r)*(1.0/beta -(1.0+delta)/b) ) - delta );
		term2 = 2.0*sqrt(omega_r)*(b-beta);
		diff = (term1+term2)/( omega + 2.0*omega_r*(1.0+beta/sqrt(omega_r)) ) ;
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


};







int main(int argc,char *argv[])
{

	double delta,delta_a,delta_rk[5], delta_a_rk[5], acc;
	double delta_i, delta_a_i;

	double a,ai,a0,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};

	int i,theory;

	FILE *fp = fopen("delta.txt","w");
	cosmo_lcdm cosmo_model_lcdm(0.3);
	cosmo_dgp cosmo_model_dgp(0.3);
	
	if(!strcmp(argv[1],"lcdm"))
	{
		printf("lcdm\n");
		theory = 0;
	}
	else
	{	printf("dgp\n");
		theory = 1;
	}

	da = 0.0001;
	ai = 0.001;
	a0 = 1.0;

	delta_i = 0.001;
	delta_a_i = .000;
	
	delta = delta_i;
	delta_a_i = delta_a_i;





	fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);

	for(a=ai;a<=a0;a+=da)
	{
		ak = a;
		delta_rk[0] = delta;
		delta_a_rk[0] = delta_a; 


		for(i=1;i<=4;++i)
		{
			if(theory==0)			
			acc  = cosmo_model_lcdm.delta_aa( ak,  delta_rk[0], delta_a_rk[0]);
			if(theory==1)			
			acc  = cosmo_model_dgp.delta_aa( ak,  delta_rk[0], delta_a_rk[0]);

			delta_a_rk[i] = da*acc;
			delta_rk[i] = da*delta_a_rk[0];

			delta_rk[0] = delta + rk_coef[i-1]*delta_rk[i];
			delta_a_rk[0] = delta_a + rk_coef[i-1]*delta_a_rk[i];

			ak = a + rk_coef[i-1]*da;


		}

		delta = delta + (1.0/6.0)*(delta_rk[1]+2.0*delta_rk[2]+2.0*delta_rk[3]+delta_rk[4]);
		delta_a = delta_a + (1.0/6.0)*(delta_a_rk[1]+2.0*delta_a_rk[2]+2.0*delta_a_rk[3]+delta_a_rk[4]);

		fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);


	}


	

}
