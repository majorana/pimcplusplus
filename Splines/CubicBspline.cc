#include "CubicBspline.h"
#include "BsplineHelper.h"
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
  P.resize(Range(0,N+2));
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
  P(N) = d(N-1)/gamma(N-1);
 
  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    P(row+1) = d(row) - mu(row)*P(row+2) - gamma(row)*P(N);
}





void
CubicBspline::Set(double start, double end, Array<double,1> &data,
		  bool interpolating, bool periodic)
{
  Periodic = periodic;
  Interpolating = interpolating;
  GridStart = start;
  GridEnd   = end;
  int N = data.size();
  
  if (!periodic) {
    GridDelta = (end-start)/(double)(data.size()-1);
    GridDeltaInv = 1.0/GridDelta;
    P.resize(data.size()+2);
    Array<double,1> d(P.size());
    d(Range(1,N)) = data;
    d(0) = 0.0;
    d(N+1) = 0.0;
    if (interpolating)
      SolveFirstDerivInterp1D (d, P);
    else {
      cerr << "Don't know how to do noninterpolating nonperiodic.\n";
      abort();
    }
  }
  else {
    GridDelta = (end-start)/(double)(data.size());
    GridDeltaInv = 1.0/GridDelta;
    P.resize(data.size()+3);

    if (interpolating) 
      SolvePeriodicInterp (data);
    else 
      P(Range(1,N)) = data;
    // Finally, assign periodic elements
    P(0)   = P(N);
    P(N+1) = P(1);
    P(N+2) = P(2);
  }
}


CubicBspline::CubicBspline()
{
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;

  dA(0,0)= 0.0; dA(0,1)= 0.0; dA(0,2)= 0.0; dA(0,3)= 0.0;
  dA(1,0)=-0.5; dA(1,1)= 1.5; dA(1,2)=-1.5; dA(1,3)= 0.5;
  dA(2,0)= 1.0; dA(2,1)=-2.0; dA(2,2)= 1.0; dA(2,3)= 0.0;
  dA(3,0)=-0.5; dA(3,1)= 0.0; dA(3,2)= 0.5; dA(3,3)= 0.0;

  d2A(0,0)= 0.0; d2A(0,1)= 0.0; d2A(0,2)= 0.0; d2A(0,3)= 0.0;
  d2A(1,0)= 0.0; d2A(1,1)= 0.0; d2A(1,2)= 0.0; d2A(1,3)= 0.0;
  d2A(2,0)=-1.0; d2A(2,1)= 3.0; d2A(2,2)=-3.0; d2A(2,3)= 1.0;
  d2A(3,0)= 1.0; d2A(3,1)=-2.0; d2A(3,2)= 1.0; d2A(3,3)= 0.0;

  d3A(0,0)= 0.0; d3A(0,1)= 0.0; d3A(0,2)= 0.0; d3A(0,3)= 0.0;
  d3A(1,0)= 0.0; d3A(1,1)= 0.0; d3A(1,2)= 0.0; d3A(1,3)= 0.0;
  d3A(2,0)= 0.0; d3A(2,1)= 0.0; d3A(1,2)= 2.0; d3A(2,3)= 0.0;
  d3A(3,0)=-1.0; d3A(3,1)= 3.0; d3A(3,2)=-3.0; d3A(3,3)= 1.0;
}
