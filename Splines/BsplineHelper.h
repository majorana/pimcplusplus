#ifndef BSPLINE_HELPER_H
#define BSPLINE_HELPER_H

#include <blitz/array.h>

using namespace blitz;

typedef enum {PERIODIC, FIXED_FIRST, FIXED_SECOND, FLAT, NATURAL} BCType;

template<typename T>
class BoundaryCondition
{
private:
  BCType BC;
  T StartDeriv, EndDeriv;
  T StartDeriv2, EndDeriv2;
public:
  inline BCType GetType() { return BC; }
  void GetDerivs (T & startDeriv, T &endDeriv) {
    assert (BC == FIXED_FIRST);
    startDeriv = StartDeriv;
    endDeriv = EndDeriv; 
  }
  void GetDerivs2 (T & startDeriv2, T &endDeriv2) {
    assert (BC == FIXED_SECOND);
    startDeriv2 = StartDeriv2;
    endDeriv2 = EndDeriv2; 
  }
  BoundaryCondition (BCType bctype)
  {
    if (bctype == FLAT) {
      BC = FIXED_FIRST;
      StartDeriv = EndDeriv = T();
    }
    else if (bctype == NATURAL) {
      BC = FIXED_SECOND;
      StartDeriv2 = EndDeriv2 = T();
    }
    else if (bctype == PERIODIC) {
      BC = PERIODIC;
    }
    else if (bctype == FIXED_FIRST) {
      cerr << "You must specify the start and end derivatives for FIXED_FIRST BC.\n";
      abort();
    }
    else if (bctype == FIXED_SECOND) {
      cerr << "You must specify the start and end derivatives for FIXED_FIRST BC.\n";
      abort();
    }
    else {
      cerr << "Unknown boundary conditions.\n";
      abort();
    }
  }
  BoundaryCondition (BCType bctype, T start, T end)
  {
    if (bctype == FIXED_FIRST) {
      StartDeriv = start;
      EndDeriv   = endl;
    }
    else if (bctype == FIXED_SECOND) {
      StartDeriv2 = start;
      EndDeriv2 = end;
    }
    else {
      cerr << "Cannot fix derivatives with periodic boundary condition.\n";
      abort();
    }
  }

};

// Solve for fixed first derivative.  
// data(0)  = dy_0/dx dx.  
// data(1)   = y_0 and
// data(M-1) = y_{M-1}.  
// data(M)   = dy_{M-1}/dx dx.  
template<typename T> inline void
SolveFirstDerivInterp1D (Array<T,1> &data, Array<T,1> &p)
{
  double ratio = 0.25;
  double ratio2 = 0.75;
  int M = data.size()-2;

  Array<double,1> d(M+2), mu(M+2);
  d = 1.5*data;
  mu = ratio;
  // First, eliminate leading coefficients
  mu(0) = 0.0;
  mu(M+1) = 0.0;
  mu(1) = ratio+ratio;
  d(0) /= (-ratio2);
  for (int row=1; row <=M; row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    mu(row) *= diagInv;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
  }
  d(M+1) /= -ratio2;
  mu(M+1) /= -ratio2;
  
  d(M+1)  -= d(M-1);
  mu(M+1) -= mu(M-1);
  d(M+1)  -= mu(M+1)*d(M);
  double diag = -1.0 - mu(M+1)*mu(M);
  p(M+1) = d(M+1)/diag;
 
  // Now go back upward, back substituting
  for (int row=M; row>=1; row--) 
    p(row) = d(row) - mu(row)*p(row+1);

  // And do 0th row
  p(0) = d(0) + p(2)*d(2);
}

#endif
