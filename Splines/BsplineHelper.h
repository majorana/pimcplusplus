#ifndef BSPLINE_HELPER_H
#define BSPLINE_HELPER_H

#include <blitz/array.h>
#include <blitz/tinymat.h>

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


// template<typename T> inline void
// SolvePeriodicInterp(Array<T,1> &data, Array<T,1> &p)
// {
//   double ratio = 0.25;
//   int N = data.size();

//   Array<double,1> d(N), gamma(N), mu(N);
//   d = 1.5*data;
//   p.resize(Range(0,N+2));
//   // First, eliminate leading coefficients
//   gamma (0) = ratio;
//   mu(0) = ratio;
//   mu(N-1) = ratio;
//   gamma(N-1) = 1.0;
//   for (int row=1; row <(N-1); row++) {
//     double diag = 1.0- mu(row-1)*ratio;
//     double diagInv = 1.0/diag;
//     gamma(row) = -ratio*gamma(row-1)*diagInv;
//     mu(row) = diagInv*ratio;
//     d(row)  = diagInv*(d(row)-ratio*d(row-1));
//     // Last row
//     d(N-1) -= mu(N-1) * d(row-1);
//     gamma(N-1) -= mu(N-1)*gamma(row-1);
//     mu(N-1) = -mu(N-1)*mu(row-1);
//   }
//   // Last row:  gamma(N-1) hold diagonal element
//   mu(N-1) += ratio;
//   gamma(N-1) -= mu(N-1)*(mu(N-2)+gamma(N-2));
//   d(N-1) -= mu(N-1) * d(N-2);
//   p(N) = d(N-1)/gamma(N-1);
 
//   // Now go back upward, back substituting
//   for (int row=N-2; row>=0; row--) 
//     p(row+1) = d(row) - mu(row)*P(row+2) - gamma(row)*p(N);
// }



// 
template<typename T> inline void
SolvePeriodicInterp1D (Array<T,1> data, Array<T,1> p)
{
  double ratio = 0.25;
  int N = data.size();

  assert (p.size() == (N+3));

  Array<double,1> d(N), gamma(N), mu(N);
  d = 1.5*data;
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
  p(N) = d(N-1)/gamma(N-1);
 
  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    p(row+1) = d(row) - mu(row)*p(row+2) - gamma(row)*p(N);
  p(0) = p(N);
  p(N+1) = p(1);
  p(N+2) = p(2);
}

inline void
SolvePeriodicInterp1D (Array<complex<double>,1> data, 
		       Array<complex<double>,1> p)
{
  int N = data.size();
  Array<double,1> dataReal(N), dataImag(N), pReal(N+3), pImag(N+3);
  for (int i=0; i<N; i++) {
    dataReal(i) = data(i).real();
    dataImag(i) = data(i).imag();
  }
  SolvePeriodicInterp1D(dataReal, pReal);
  SolvePeriodicInterp1D(dataImag, pImag);

  for (int i=0; i<N+3; i++) 
    p(i) = complex<double>(pReal(i), pImag(i));
}



// We multiply both sides by 1/4 to make diagonals equal 1
template<typename T> inline void
SolveDerivInterp1D (Array<T,1> data, Array<T,1> p,
		    TinyVector<double,4> abcdInitial,
		    TinyVector<double,4> abcdFinal)
{
  assert (p.size() == (data.size()+2));
  double al = 0.25*abcdInitial[0];
  double bl = 0.25*abcdInitial[1];
  double cl = 0.25*abcdInitial[2];
  double dl = 1.5 *abcdInitial[3];
  double ar = 0.25*abcdFinal[0];
  double br = 0.25*abcdFinal[1];
  double cr = 0.25*abcdFinal[2];
  double dr = 1.5 *abcdFinal[3];
    
  double ratio = 0.25;
  int M = data.size();

  Array<double,1> d(M+2), mu(M+2);
  d(Range(1,M)) = 1.5*data;
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


inline void
SolveDerivInterp1D (Array<complex<double>,1> data, 
		    Array<complex<double>,1> p,
		    TinyVector<double,4> abcdInitial,
		    TinyVector<double,4> abcdFinal)
{
  int N = data.size();
  Array<double,1> dataReal(N), dataImag(N), pReal(N+3), pImag(N+3);
  for (int i=0; i<N; i++) {
    dataReal(i) = data(i).real();
    dataImag(i) = data(i).imag();
  }
  SolveDerivInterp1D(dataReal, pReal, abcdInitial, abcdFinal);
  SolveDerivInterp1D(dataImag, pImag, abcdInitial, abcdFinal);

  for (int i=0; i<N+3; i++) 
    p(i) = complex<double>(pReal(i), pImag(i));
}


template <typename GridType>
class NUBsplineBasis
{
private:
  GridType &grid;
  inline void BuildBasis();
  Array<double,1> xVals;
  Array<TinyVector<double,3>,1> dxInv;
public:
  inline void operator() (double x, TinyVector<double,4>& bfuncs);
  inline void Evaluate   (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv);
  inline void Evaluate   (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv,
			  TinyVector<double,4> &deriv2);
  NUBsplineBasis (GridType &newGrid) : grid(newGrid)
  {
    BuildBasis();
  }
  
};

template<typename t> void
NUBsplineBasis<t>::BuildBasis()
{
  int N = grid.NumPoints;
  dxInv.resize(grid.NumPoints+2);
  xVals.resize(N+4);
  for (int i=0; i<N; i++)
    xVals(i+2) = grid(i);

  xVals(0) = grid(0) - 2.0*(grid(1)-grid(0));
  xVals(1) = grid(0) - 1.0*(grid(1)-grid(0));
  xVals(N+2) = grid(N-1) + 1.0*(grid(N-1)-grid(N-2));
  xVals(N+3) = grid(N-1) + 2.0*(grid(N-1)-grid(N-2));
  for (int i=0; i<N+2; i++) {
    for (int j=0; j<3; j++) {
      if ((i-2) <= 0.0)
	dxInv(i)[j] = 1.0/(grid(j+1)-grid(0));
      else if (i+j-1 >= N)
	dxInv(i)[j] = 1.0/(grid(N-1)-grid(N-j-2));
      else //((i-2)>=0) && (i+j-1)<grid.NumPoints)
	dxInv(i)[j] = 1.0/(grid(i+j-1)-grid(i-2));
//       else
// 	dxInv(i)[j] = 0.0;
    }
  }
}


template<typename t> void
NUBsplineBasis<t>::operator()(double x, TinyVector<double,4> &bfuncs)
{
  double b1[2], b2[3];
  int i = grid.ReverseMap (x);
  int i2 = i+2;
//   b1[0] = (grid(i+1)-x)/(grid(i+1)-grid(i));
//   b1[1] = (x-grid(i))/(grid(i+1)-grid(i));
//   b2[0] = (grid(i+1)-x)/(grid(i+1)-grid(i-1))*b1[0];
//   b2[1] = ((x-grid(i-1))/(grid(i+1)-grid(i-1))*b1[0] +
// 	   (grid(i+2)-x)/(grid(i+2)-grid(i))*b1[1]);
//   b2[2] = (x-grid(i))/(grid(i+2)-grid(i))*b1[1];
//   bfuncs[0] = (grid(i+1)-x)/(grid(i+1)-grid(i-2)) * b2[0];
//   bfuncs[1] = ((x-grid(i-2))/(grid(i+1)-grid(i-2))*b2[0] +
// 	       (grid(i+2)-x)/(grid(i+2)-grid(i-1))*b2[1]);
//   bfuncs[2] = ((x-grid(i-1))/(grid(i+2)-grid(i-1))*b2[1] +
// 	       (grid(i+3)-x)/(grid(i+3)-grid(i))*b2[2]);
//   bfuncs[3] = (x-grid(i))/(grid(i+3)-grid(i))*b2[2];

//   b1[0]     = (grid(i+1)-x)  * dxInv(i+2)[0];
//   b1[1]     = (x-grid(i))    * dxInv(i+2)[0];
//   b2[0]     = (grid(i+1)-x)  * dxInv(i+1)[1] * b1[0];
//   b2[1]     = ((x-grid(i-1)) * dxInv(i+1)[1] * b1[0]+
// 	       (grid(i+2)-x) * dxInv(i+2)[1] * b1[1]);
//   b2[2]     = (x-grid(i))    * dxInv(i+2)[1] * b1[1];
//   bfuncs[0] = (grid(i+1)-x)  * dxInv(i  )[2] * b2[0];
//   bfuncs[1] = ((x-grid(i-2)) * dxInv(i  )[2] * b2[0] +
// 	       (grid(i+2)-x) * dxInv(i+1)[2] * b2[1]);
//   bfuncs[2] = ((x-grid(i-1)) * dxInv(i+1)[2] * b2[1] +
// 	       (grid(i+3)-x) * dxInv(i+2)[2] * b2[2]);
//   bfuncs[3] = (x-grid(i))    * dxInv(i+2)[2] * b2[2];

  b1[0]     = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]     = (x-xVals(i2))    * dxInv(i+2)[0];
  b2[0]     = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]     = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
	       (xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]     = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  bfuncs[0] = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
	       (xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
	       (xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3] = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];
}




#endif
