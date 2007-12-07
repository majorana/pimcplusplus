#ifndef VINET_EOS2_H
#define VINET_EOS2_H

#include <blitz/tinyvec.h>

using namespace blitz;

class VinetEOSClass
{
private:
  double V0, E0, B0, B0p;
  double third;
  double Atomic2GPa;
public:
  inline double operator()(double V);
  inline TinyVector<double,4> Grad(double V);
  inline TinyVector<double,4> GradFD(double V);
  inline double Pressure(double V);
  inline double PressureFD(double V);
  inline void SetParams (TinyVector<double,4> params);
  inline TinyVector<double,4> GetParams ();
  inline double GetB0p();
  VinetEOSClass() : third(1.0/3.0), Atomic2GPa(29421.01)
  { 
  }
};

inline void
VinetEOSClass::SetParams (TinyVector<double,4> p)
{
  V0  = p[0];
  E0  = p[1];
  B0  = p[2];
  B0p = p[3];
}

inline TinyVector<double,4>
VinetEOSClass::GetParams ()
{
  TinyVector<double,4> p;
  p[0] = V0;
  p[1] = E0;
  p[2] = B0;
  p[3] = B0p;
  return p;
}



inline double
VinetEOSClass::operator()(double V)
{
  double x = cbrt(V/V0);
  double eta = 1.5*(B0p-1.0);
  
  return -2.0*B0*exp(eta*(1.0-x))/(x*x*(B0p-1.0)*(B0p-1.0))*(3.0*(B0p-1.0)*V + (5.0-3.0*B0p)*x*x*V0) + E0;
}


inline double
VinetEOSClass::Pressure(double V)
{
  double x = cbrt(V/V0);
  double eta = 1.5*(B0p-1.0);
  double dE_dV = -3.0*B0/(x*x)*(1.0-x)*exp(eta*(1.0-x));
  return -Atomic2GPa * dE_dV;
}


inline double
VinetEOSClass::PressureFD(double V)
{
  return Atomic2GPa*((*this)(V+1.0e-5)-(*this)(V-1.0e-5))/-2.0e-5;
}

inline double
VinetEOSClass::GetB0p()
{
  return B0p;
}


inline TinyVector<double,4>
VinetEOSClass::VinetEOSClass::Grad (double V)
{
  double x = cbrt(V/V0);
  double eta = 1.5*(B0p-1.0);

  double expval = exp(eta*(1.0-x));
  double dE_dE0 = 1.0;
  double dE_dV0 = -B0/(V0*x*x*(B0p-1.0)*(B0p-1.0))*expval*
    (3.0*(B0p-1.0)*V*(3.0+B0p*(x-1.0) - x) + 2.0*(5.0-3.0*B0p)*x*x*V0);
  double dE_dB0 = -2.0*expval/(x*x)*(3.0*(B0p-1.0)*V + (5.0-3.0*B0p)*x*x*V0)/
    ((B0p-1.0)*(B0p-1.0));
  double dE_dB0p = B0/((B0p-1.0)*(B0p-1.0)*(B0p-1.0))*expval/(x*x)*
    (3.0*(B0p-1.0)*V*(10.0+3.0*B0p*(x-2.0) - 3.0*x)+(29.0-30.0*B0p + 9.0*B0p*B0p)
     *x*x*V0);
  
  return TinyVector<double,4> (dE_dV0, dE_dE0, dE_dB0, dE_dB0p);
}

inline TinyVector<double,4>
VinetEOSClass::GradFD(double V)
{
  TinyVector<double,4> p, p_plus, p_minus, grad;
  double Eplus, Eminus;
  p = GetParams();
  for (int i=0; i<4; i++) {
    p_plus = p; p_minus = p;
    p_plus[i]  += 1.0e-7;
    p_minus[i] -= 1.0e-7;
    SetParams(p_plus);   Eplus  = (*this)(V);
    SetParams(p_minus);  Eminus = (*this)(V);
    grad[i] = (Eplus-Eminus)/(2.0e-7);
  }
  return grad;
}




#endif
