#ifndef OPTIMIZED_BREAKUP_H
#define OPTIMIZED_BREAKUP_H

#include <blitz/array.h>

using namespace blitz;

class BasisClass
{
protected:
  double r_c;
public:
  /// Set the cutoff radius
  virtual void Set_rc(double rc) = 0;
  inline double Get_rc() { return r_c; }
  /// Returns the number of basis elements
  virtual int NumElements() = 0;
  /// Returns the basis element n evaluated in real space at r
  virtual double h(int n, double r) = 0;
  /// Returns the basis element n evaluated in k space at k
  virtual double c(int n, double k) = 0;
  BasisClass() : r_c (0.0)
  { /* do nothing */ }
};



/// Locally Piecewise Quintic Hermite Interpolant
class LPQHI_BasisClass : public BasisClass
{
private:
  int NumKnots;
  double delta, deltaInv;
  TinyMatrix<double,3,6> S;
public:
  // n must be at least 2;
  void SetNumKnots(int n);
  void Set_rc(double rc);
  double Get_rc();
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


#endif
