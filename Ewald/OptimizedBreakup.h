#ifndef OPTIMIZED_BREAKUP_H
#define OPTIMIZED_BREAKUP_H

#include <blitz/array.h>

using namespace blitz;

class BasisClass
{
public:
  //protected:
  double r_c;
  TinyVector<double,3> Box;
  double Omega;
public:
  /// Set the cutoff radius
  virtual void Set_rc(double rc) = 0;
  inline double Get_rc() { return r_c; }
  inline void SetBox (TinyVector<double,3> box);
  inline TinyVector<double,3> GetBox ();
  /// Returns the number of basis elements
  virtual int NumElements() = 0;
  /// Returns the basis element n evaluated in real space at r
  virtual double h(int n, double r) = 0;
  /// Returns the basis element n evaluated in k space at k
  virtual double c(int n, double k) = 0;
  double c_numerical (int n, double k);
  /// This returns the coefficent of the nth basis function
  //virtual double  Get_t(int n) const     = 0;
  /// This sets the coefficent of the nth basis function
  //virtual void    Set_t(int n, double t) = 0;
  /// This returns the linear combination of the basis functions with
  /// coefficients t_n
  //virtual double f (double r) = 0;
  BasisClass() : r_c (0.0)
  { /* do nothing */ }
};


class OptimizedBreakup
{
private:
  BasisClass &Basis;
  void Addk(double k);
public:
  // First element is |k|, second is degeneracy of the point.
  Array<TinyVector<double,2>,1> kpoints;
  void SetkVecs(double kc, double kcont, double kMax);
  /// kc is the k-space cutoff for the Ewald sum.  
  /// kMax is largest k we use in determining the error in the breakup.  
  /// t is the set of coefficients of the breakup.
  /// inFit is a boolean array telling whether t_n should be optimized
  /// or left at its initial value.
  void DoBreakup (const Array<double,1> &Vk, Array<double,1> &t, 
		  const Array<bool,1> &adjust);
  /// Same as above, but we assume that all t's are adjusted.
  void DoBreakup (const Array<double,1> &Vk, Array<double,1> &t);
  OptimizedBreakup (BasisClass &basis) : Basis(basis)
  { /* Do nothing */ }
};


inline void BasisClass::SetBox (TinyVector<double,3> box)
{
  Box = box;
  Omega = box[0]*box[1]*box[2];
}

inline TinyVector<double,3> BasisClass::GetBox ()
{
  return Box;
}



/// Locally Piecewise Quintic Hermite Interpolant
class LPQHI_BasisClass : public BasisClass
{
  /// public is HACK
  //private:
public:
  int NumKnots;
  double delta, deltaInv;
  TinyMatrix<double,3,6> S;
  /// The following are helpers to calculate the Fourier tranform of
  /// the basis functions
  inline complex<double> Eplus(int i, double k, int n);
  inline complex<double> Eminus(int i, double k, int n);
  inline double Dplus(int i, double k, int n);
  inline double Dminus(int i, double k, int n);
  Array<double,1> tvec;
public:
  // n must be at least 2;
  void SetNumKnots(int n);
  void Set_rc(double rc);
  int NumElements();
  double h(int n, double r);
  double c(int n, double k);
  LPQHI_BasisClass() : NumKnots(0), delta(0.0) 
  { 
    S = 
      1.0,   0.0,   0.0, -10.0,  15.0,  -6.0,
      0.0,   1.0,   0.0,  -6.0,   8.0,  -3.0,
      0.0,   0.0,   0.5,  -1.5,   1.5,  -0.5;
  }
};

inline complex<double> LPQHI_BasisClass::Eplus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0);

  if (n == 0) {
    complex<double> e1(cos(k*delta)-1.0, sin(k*delta));
    complex<double> e2(cos(k*delta*i),   sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else {
    complex<double> t1, t2;
    double sign = 1.0;
    t1 = complex<double>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    t2=-(double)n/delta*Eplus(i,k,n-1);;
    return (-(eye/k)*(t1+t2));
  }
}

inline complex<double> LPQHI_BasisClass::Eminus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0);

  if (n == 0) {
    complex<double> e1(cos(k*delta)-1.0, -sin(k*delta));
    complex<double> e2(cos(k*delta*i),    sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else {
    complex<double> t1, t2;
    double sign = (n & 1) ? -1.0 : 1.0;
    t1 = sign*
      complex<double> (cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    t2=-(double)n/delta*Eminus(i,k,n-1);
    return (-(eye/k)*(t1+t2));
  }
}


inline double LPQHI_BasisClass::Dplus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0); 
  complex<double> Z1 = Eplus(i,k,n+1);
  complex<double> Z2 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag());
}


inline double LPQHI_BasisClass::Dminus(int i, double k, int n)
{
  complex<double> eye(0.0, 1.0); 
  complex<double> Z1 = Eminus(i,k,n+1);
  complex<double> Z2 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag());
}

#endif
