#ifndef FREE_PARTICLES_H
#define FREE_PARTICLES_H

double LogFreeIntegrand(int l, double x, double xp, double delta, 
			double lambda, double beta);
double LogFreeIntegrand_Hermite(int l, double x, double xp, double delta, 
				double lambda, double beta);

double drho0_dbeta_l        (int l, double x, double xp, 
		             double lambda, double beta);
double drho0_dbeta_l_scaled (int l, double x, double xp, 
			     double lambda, double beta);

double drho0_dbeta_l_FD(int l, double x, double xp, double lambda,double beta);

double Log_rho0_l (int l, double x, double xp, double lambda, double beta);

double drho0_dx_l(int l, double x, double xp, double lambda, double beta);
double d2rho0_dx2_l(int l, double x, double xp, double lambda, double beta);
void CheckLogFreeIntegrand();

double Log_rho0 (double r, double rp, double costheta, 
		 double lambda, double beta);
double rho0     (double r, double rp, double costheta,
	         double lambda, double beta);

double drho0_dbeta_scaled(double r, double rp, double costheta, 
			  double lambda, double beta);
double drho0_dbeta       (double r, double rp, double costheta, 
		          double lambda, double beta);

void Test_Rho0_Derivs();

#endif
