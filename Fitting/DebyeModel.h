#ifndef DEBYE_MODEL_H
#define DEBYE_MODEL_H

#include "../Blitz.h"
#include "gsl_sf_debye.h"

class DebyeModel_T
{
private:
  double Theta, ThetaInv;
  const double kB; // In Hartrees per Kelvin
  double N;
public:
  inline void SetN    (int n);
  inline void SetTheta(double theta);
  inline double F(double T);
  void Fit (A
  DebyeModel_T() : kB(3.16681526543384e-06)
  {
  }
};

inline void
DebyeModel_T::SetN(int n)
{
  N = (double)n;
}

inline void
DebyeModel_T::SetTheta (double theta)
{
  Theta = theta; 
  ThetaInv = 1.0/theta;
}

inline double
DebyeModel_T::F(double T)
{
  double x = T * ThetaInv;
  
  double D = gsl_sf_debye_3(x);
  return 3.0*N*kB*T*log1p(-exp(-1.0/x)) - N*kBT*D;
}

#endif
