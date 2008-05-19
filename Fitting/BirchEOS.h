#ifndef BIRCH_EOS_H
#define BIRCH_EOS_H

#include <blitz/tinyvec.h>

using namespace blitz;

template<int N>
class BirchEOSClass
{
private:
  TinyVector<double,N> bn;
  double V0, E0;
  double third;
  double Atomic2GPa;
public:
  inline double operator()(double V);
  inline TinyVector<double,4> Grad(double V);
  inline TinyVector<double,4> GradFD(double V);
  inline double Pressure(double V);
  inline double PressureFD(double V);
  inline void SetParams (TinyVector<double,N+2> params);
  inline TinyVector<double,4> GetParams ();
  inline double GetB0p();
  inline double K_T (double V);

  BirchEOSClass() : third(1.0/3.0), Atomic2GPa(29421.01)
  { 
  }
};

template<int N>
inline void
BirchEOSClass<N>::SetParams (TinyVector<double,N+2> p)
{
  V0  = p[0];
  E0  = p[1];
  for (int i=0; i<N; i++)
    bn[i] = p[i+2];
}

template<int N>
inline TinyVector<double,N+2>
BirchEOSClass<N>::GetParams ()
{
  TinyVector<double,N+2> p;
  p[0] = V0;
  p[1] = E0;
  for (int i=0; i<N; i++)
    p[i+2] = bn[i];
  return p;
}


template<int N>
inline double
BirchEOSClass<N>::operator()(double V)
{
  double x = cbrt(V/V0);
  double t = x*x-1;
  double t_np2 = t*t;
  double E = E0;
  for (int n=0; n<N; n++) {
    E += bn[n] * t_np2;
    t_np2 *= t;
  }
  return E;
}


template<int N>
inline double
BirchEOSClass<N>::Pressure(double V)
{
  double x = cbrt(V/V0);
  double t = x*x-1;
  double tprime = -2.0/(3.0*V)*(t + 1.0);
  double P = 0.0;
  double t_np1 = t;
  for (int n=0; n<N; n++) {
    P += (double)(n+2)*bn[n]*t_np1;
    t_np1 *= t;
  }
  P *= 2.0/(3.0*V)*(t + 1.0);

  return Atomic2GPa * dE_dV;
}

template<int N>
inline double
BirchEOSClass<N>::K_T (double V)
{
  double eps = 1.0e-7;
  return -V*(Pressure(V+eps)-Pressure(V-eps))/(2.0*eps);
}

template<int N>
inline double
BirchEOSClass<N>::PressureFD(double V)
{
  return Atomic2GPa*((*this)(V+1.0e-5)-(*this)(V-1.0e-5))/-2.0e-5;
}


template<int N>
inline double
BirchEOSClass<N>::GetB0p()
{
}


template<int N>
inline TinyVector<double,N+2>
BirchEOSClass::BirchEOSClass<N>::Grad (double V)
{
}

template<int N>
inline TinyVector<double,N+2>
BirchEOSClass<N>::GradFD(double V)
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
