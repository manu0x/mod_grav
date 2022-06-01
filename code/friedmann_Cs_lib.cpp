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
		//printf("H0 %lf\n",H0);

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
		//printf("H0 %lf  B1  %lf\n",H0,B1);

	}

	double Hsqr(double a,double a0=1.0)
	{
		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		val = H0*H0*(   0.5*omega + B0/6.0  + sqrt((0.5*omega + B0/6.0)*(0.5*omega + B0/6.0) + B1*B1/3.0 ) );
		return(val);


	}

	double HbyH0(double a,double a0=1.0)
	{
		double val,omega;
		omega = omega_dm_0*pow(a0/a,3.0);
		val = sqrt(   0.5*omega + B0/6.0  + sqrt((0.5*omega + B0/6.0)*(0.5*omega + B0/6.0) + B1*B1/3.0 ) );
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





int f_sigma_log_like_bimetric(double omega_dm0,double B1,double sig8,double *data,double *fs8,double *ratio_den,long int *cn,int n,int ax=1,int print_file=0,double a0=1.0)
{
	double D[3],D_a[3],D_rk[5][3], D_a_rk[5][3], acc[3];
	double D_i[3], D_a_i[3];

	double a,ai,ai_burn,ak,da;
	double rk_coef[4] = {0.5,0.5,1.0,1.0};
	double om_dm_0, model_param;
	double D1_0;
	double a_lst[n];
	double H_H0[n];	
	int a_lst_cntr=0;
	double a_lst_now ;
	
	double log_like;
	

	int i,j,cn_i;

	//FILE *fp_fs8 = fopen("bimetric_fs8_r.txt","w");

	om_dm_0 = omega_dm0;

	for(i=0,j=0;i<n;++i,++j)
	{
		if(ax)
		a_lst[i] = data[j];
		
		else
		a_lst[i] = a0/(1.0+data[j]);


	}
	
	model_param = B1;
	
	cosmo_bigravity cosmo_model_bigravity(om_dm_0,model_param);


	da = 0.0001;
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
	a_lst_now = a_lst[0];

	//fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",a/ai,delta,delta_a,delta/delta_i,delta_a/delta_a_i);

	for(a=ai_burn,cntr=0;a<=a0;a+=da,++cntr)
	{
		ak = a;

		if((a>=a_lst_now)&&(a_lst_cntr<n))
		{
			
			for(cn_i=0;cn_i<cn[a_lst_cntr];++cn_i)
			{fs8[a_lst_cntr+cn_i] = sig8*a*D_a[0];

			 H_H0[a_lst_cntr+cn_i] = cosmo_model_bigravity.HbyH0( a_lst[a_lst_cntr+cn_i]);
			}
			
			//fprintf(fp_fs8,"%d\t%d\t%lf\n",a_lst_cntr,cn[a_lst_cntr],fs8[a_lst_cntr]);	

			a_lst_cntr+=cn[a_lst_cntr];
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
		
			
		}



		D_rk[0][0] = D[0];
		D_a_rk[0][0] = D_a[0]; 

		D_rk[0][1] = D[1];
		D_a_rk[0][1] = D_a[1]; 

		D_rk[0][2] = D[2];
		D_a_rk[0][2] = D_a[2]; 
		
		

		for(i=1;i<=4;++i)
		{
			
			cosmo_model_bigravity.pert_delta_aa(acc, D_rk[0], D_a_rk[0],a);

			
		   for(j=0;j<3;++j)
			{D_a_rk[i][j] = da*acc[j];
			 D_rk[i][j] = da*D_a_rk[0][j];
				

			 

			 D_rk[0][j] = D[j] + rk_coef[i-1]*D_rk[i][j];
			 D_a_rk[0][j] = D_a[j] + rk_coef[i-1]*D_a_rk[i][j];
			}

			ak = a + rk_coef[i-1]*da;
			
		}


	    for(j=0;j<3;++j)
		{
		 D[j] = D[j] + (1.0/6.0)*(D_rk[1][j]+2.0*D_rk[2][j]+2.0*D_rk[3][j]+D_rk[4][j]);
		 D_a[j] = D_a[j] + (1.0/6.0)*(D_a_rk[1][j]+2.0*D_a_rk[2][j]+2.0*D_a_rk[3][j]+D_a_rk[4][j]);

		

		}


	}

	D1_0 = D[0];

	for(i=0,j=0;i<n;++i,++j)
	{

		fs8[i] = (fs8[i]/D1_0)*(H_H0[i]/ratio_den[i]);
	/*
		if(ax)
		{if(print_file)
		 fprintf(fp_fs8,"%lf\t%lf\t%lf\n",a_lst[i],fs8[i],(H_H0[i]/ratio_den[i]));
		}
		else
		{if(print_file)
		  fprintf(fp_fs8,"%lf\t%lf\t%lf\n",data[j],fs8[i],(H_H0[i]/ratio_den[i]));

		}
	*/
	}


	//fclose(fp_fs8);
	
	return(1);


}


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



extern "C" {
    int bimetric_fs8_log_like(double omega_dm0,double B1,double sig8,double *data,double *fs8,double *ratio_den,long int *cn,int n,int ax=1,int print_file=0,double a0=1.0)
    {
        return f_sigma_log_like_bimetric(omega_dm0,B1,sig8,data,fs8,ratio_den,cn, n, ax, print_file, a0);

	
    }
}











///##########################################################	
