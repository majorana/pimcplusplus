#include "CubicBspline.h"
#include <iostream>
#include <cstdlib>

// Solve tridiagonal linear system in periodic boundary conditions
void
CubicBspline::SolvePeriodicInterp(Array<double,1> &data)
{
  double ratio = 0.25;
  int N = data.size();

  Array<double,1> d(N), gamma(N), mu(N);
  d = 1.5*data;
  P.resize(N);
  // First, eliminate leading coefficients
  gamma (0) = ratio;
  mu(0) = ratio;
  mu(N-1) = ratio;
  gamma(N-1) = 1.0;
  for (int row=1; row <(N-1); row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    gamma(row) = -ratio*gamma(row-1)*diagInv;
    mu(row) = diagInv*ratio;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
    // Last row
    d(N-1) -= mu(N-1) * d(row-1);
    gamma(N-1) -= mu(N-1)*gamma(row-1);
    mu(N-1) = -mu(N-1)*mu(row-1);
  }
  // Last row:  gamma(N-1) hold diagonal element
  mu(N-1) += ratio;
  gamma(N-1) -= mu(N-1)*(mu(N-2)+gamma(N-2));
  d(N-1) -= mu(N-1) * d(N-2);
  P(N-1) = d(N-1)/gamma(N-1);
 
  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    P(row) = d(row) - mu(row)*P(row+1) - gamma(row)*P(N-1);
}

void
CubicBspline::Set(double start, double end, Array<double,1> &data,
		  bool interpolating, bool periodic)
{
  Periodic = periodic;
  Interpolating = interpolating;
  GridStart = start;
  GridEnd   = end;
  GridDelta = (end-start)/(double)(data.size());
  GridDeltaInv = 1.0/GridDelta;
  P.resize(data.size());
  
  if (!periodic) {
    cerr << "Non-periodic B-splines not yet implemented.\n";
    abort();
  }
  else {
    if (interpolating) 
      SolvePeriodicInterp (data);
    else 
      P = data;
  }
}


CubicBspline::CubicBspline()
{
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;
}
