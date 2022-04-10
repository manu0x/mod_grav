using namespace std;

#include <math.h>
#include <stdio.h>

class cosmo
{

	private:
	double omega_dm_0, H0;
	int model;

	public:
	cosmo(omega_dm_0_val=0.3,h_val=0.7,model=0)
	{
		omega_dm_0 = omega_dm_0_val;
		H0 = h_val;
		model = model;

	}


	double Hb_t_by_Hb2(double a)
	{

		return(1.0);



	}

	double H_Diff(double a, double delta)
	{

	return(1.0);

	}

	double delta_aa(double a, double delta, double delta_a)
	{
		double HbtbyHb2, diff;

		HbtbyHb2 = Hb_t_by_Hb2(a);
		diff = H_Diff(a, delta);

		acc = (1.0/a)*(3.0 + HbtbyHb2)*delta_da + (4.0/3.0)*delta_a*delta_a/(1.0+delta) - 3.0*(1.0+delta)*diff;

	}


};










int main()
{

	double delta,delta_a,delta_rk[5], delta_a_rk[5], acc;

	double a,ai,a0,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};

	int i;
	
	cosmo cosmo_model();

	for(a=ai;a<=a0;a+=da)
	{
		ak = a;
		delta_rk[0] = delta;
		delta_a_rk[0] = delta_a; 

		for(i=1;i<=4;++i)
		{
			acc  = cosmo_model.delta_aa( ak,  delta_rk[0], delta_a_rk[0]);

			delta_a_rk[i] = da*acc;
			delta_rk[i] = da*delta_a_rk[0];

			delta_rk[0] = delta + rk_coef[i-1]*delta_rk[i];
			delta_a_rk[0] = delta_a + rk_coef[i-1]*delta_a_rk[i];

			ak = a + rk_coef[i-1]*da;


		}

		delta = delta + (1.0/6.0)*(delta_rk[1]+2.0*delta_rk[2]+2.0*delta_rk[3]+delta_rk[4]);
		delta_a = delta_a + (1.0/6.0)*(delta_a_rk[1]+2.0*delta_a_rk[2]+2.0*delta_a_rk[3]+delta_a_rk[4]);


	}




}
