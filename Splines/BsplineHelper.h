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
  T Val;
public:
  inline BCType GetType() { return BC; }
  T GetVal () {
    assert (BC == FIXED_FIRST || BC == FIXED_SECOND);
    return Val;
  }
  BoundaryCondition (BCType bctype)
  {
    if (bctype == FLAT) {
      BC = FIXED_FIRST;
      Val = T();
    }
    else if (bctype == NATURAL) {
      BC = FIXED_SECOND;
      Val = T();
    }
    else if (bctype == PERIODIC) {
      BC = PERIODIC;
    }
    else if (bctype == FIXED_FIRST) {
      cerr << "You must specify the derivative value for FIXED_FIRST BC.\n";
      abort();
    }
    else if (bctype == FIXED_SECOND) {
      cerr << "You must specify the FIXED_SECOND BC.\n";
      abort();
    }
    else {
      cerr << "Unknown boundary conditions.\n";
      abort();
    }
  }
  BoundaryCondition (BCType bctype, T val)
  {
    if ((bctype != FIXED_FIRST) && (bctype != FIXED_SECOND)) {
      cerr << "Cannot fix derivatives with periodic, flat, or natural boundary conditions.\n";
      abort();
    }
    Val = val;
    BC = bctype;
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
  p(0) = d(0) + p(2)/**d(2)*/;
}

// We multiply both sides by 1/4 to make diagonals equal 1
template<typename T> inline void
SolveDerivInterp1D (Array<T,1> &data, Array<T,1> &p,
		    TinyVector<T,4> abcdInitial,
		    TinyVector<T,4> abcdFinal)
{
  double al = 0.25*abcdInitial[0];
  double bl = 0.25*abcdInitial[1];
  double cl = 0.25*abcdInitial[2];
  double dl = 1.5 *abcdInitial[3];
  double ar = 0.25*abcdFinal[0];
  double br = 0.25*abcdFinal[1];
  double cr = 0.25*abcdFinal[2];
  double dr = 1.5 *abcdFinal[3];
    
  double ratio = 0.25;
  int M = data.size()-2;

  Array<double,1> d(M+2), mu(M+2);
  d = 1.5*data;
  mu = ratio;
  // First, eliminate leading coefficients
  double alInv = 1.0/al;
  bl *= alInv;
  cl *= alInv;
  dl *= alInv;

  d(0) = dl;  
  mu(0) = bl;
  mu(1) = ratio - ratio*cl;

  for (int row=1; row <=M; row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    mu(row) *= diagInv;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
  }
  
  br -= ar*mu(M-1);
  dr -= ar*d(M-1);
  cr -= br*mu(M);
  dr -= br*d(M);
  p(M+1) = dr/cr;
   
  // Now go back upward, back substituting
  for (int row=M; row>=1; row--) 
    p(row) = d(row) - mu(row)*p(row+1);

  // And do 0th row
  p(0) = dl -bl*p(1) - cl*p(2);
}

#endif
