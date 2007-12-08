#ifndef DEBYE_MODEL_H
#define DEBYE_MODEL_H

#include "../Blitz.h"
#include "gsl/gsl_sf_debye.h"
#include <vector>

class DebyeModel
{
private:
  double Theta, ThetaInv;
  double kB; // In Hartrees per Kelvin
  double N;
  double Chi2 (double theta, Array<double,1> &Fvals,
	       Array<double,1> &Tvals);
  double MinTheta (Array<double,1> &Fvals, Array<double,1> &Tvals,
		   double startTheta, double endTheta);
public:
  inline void SetN    (int n);
  inline void SetTheta(double theta);
  inline double F(double T);
  inline double dF_dTheta (double T);
  double OptTheta (Array<double,1> &Fvals, Array<double,1> &Tvals);
  inline double CalcTheta (double F, double T);
  DebyeModel() : kB(3.16681526543384e-06)
  {
  }
};


class DebyeFreeEnergy
{
private:
  vector<double> V, Theta, F0;
  DebyeModel Debye;
  Array<double,1> ThetaCoefs;
public:
  void AddModel (double V, vector<double> &F,
		 vector<double> &T);
  inline double Theta_V (double V);
  inline double dTheta_dV (double V);
  
  void FitTheta_V(int numCoefs);
  // Helmholtz free energy
  double F(double V, double T);
  // Thermal pressure
  double P(double V, double T);
  double P_FD(double V, double T);

  DebyeFreeEnergy()
  {
    Debye.SetN(2);
  }
};

inline double 
DebyeFreeEnergy::Theta_V(double V)
{
  double V2i = 1.0;
  double tsum = 0.0;
  for (int i=0; i<ThetaCoefs.size(); i++) {
    tsum += ThetaCoefs(i) * V2i;
    V2i *= V;
  }
  return tsum;
}

inline double 
DebyeFreeEnergy::dTheta_dV(double V)
{
  double V2im1 = 1.0;
  double tsum = 0.0;
  for (int i=1; i<ThetaCoefs.size(); i++) {
    tsum += ThetaCoefs(i) * (double)i * V2im1;
    V2im1 *= V;
  }
  return tsum;
}


inline double
DebyeModel::CalcTheta(double Fval, double T)
{
  double thetaLow = 0.0; 
  double thetaHigh = 10000.0;

  while ((thetaHigh - thetaLow) > 1.0e-10) {
    double theta = 0.5*(thetaHigh + thetaLow);
    SetTheta (theta);
    double Ftry = F(T);
    if (Ftry > Fval) 
      thetaLow = theta;
    else
      thetaHigh = theta;
  }
  return 0.5*(thetaHigh + thetaLow);
}


inline void
DebyeModel::SetN(int n)
{
  N = (double)n;
}

inline void
DebyeModel::SetTheta (double theta)
{
  Theta = theta; 
  ThetaInv = 1.0/theta;
}

inline double
DebyeModel::F(double T)
{
  double x = T * ThetaInv;
  if (x == 0.0)
    x = 1.0e-10;
  
  double D = gsl_sf_debye_3(1.0/x);
  double f = N*kB*T*(3.0*log1p(-exp(-1.0/x)) - D  /*-9.0/(8.0*x) */);
  // cerr << "theta = " << Theta << "  T = " << T << "  F = " << f << endl;
//   fprintf (stderr, "theta = %1.5f   T = %1.5f   F = %1.8e  D = %1.8e\n", 
// 	   Theta, T, F);
  return f;
}

inline double
DebyeModel::dF_dTheta (double T)
{
  double x = T * ThetaInv;
  double D3 = gsl_sf_debye_3(1.0/x);
  
  return ThetaInv * 3.0*N*kB*T*D3;
}


#endif
