#ifndef BSPLINE_HELPER_H
#define BSPLINE_HELPER_H

#include <blitz/array.h>
#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include <blitz/tinyvec-et.h>
#include <complex>
//#include "../Blitz.h"

using namespace blitz;
using namespace std;

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
  GridType *GridPtr;

  // xVals is just the grid points, augmented by two extra points on
  // either side.  These are necessary to generate enough basis
  // functions. 
  Array<double,1> xVals;
  // dxInv(i)[j] = 1.0/(grid(i+j-1)-grid(i-2))
  // This is used to avoid division in evalating the splines
  Array<TinyVector<double,3>,1> dxInv;
  bool Periodic;
public:
  // Initialized xVals and dxInv.
  inline double Grid (int i) { return ((*GridPtr)(i)); }
  inline GridType& GetGrid() { return (*GridPtr); }
  inline void Init(GridType *gridPtr, bool periodic=false);
  // Evaluates the basis functions at a give value of x.  Returns the
  // index of the first basis function
  inline int  operator() (double x, TinyVector<double,4>& bfuncs) const;
  // Same as above, but also computes first derivatives
  inline int  operator() (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv) const;
  // Same as above, but also computes second derivatives
  inline int  operator() (double x, TinyVector<double,4>& bfunc,
			  TinyVector<double,4> &deriv,
			  TinyVector<double,4> &deriv2) const;
  // These versions take the grid point index.  These are needed for
  // the boundary conditions on the interpolating equations.  The
  // complex version return the basis functions in both the real and
  // imaginary parts.
  inline void operator() (int i, TinyVector<double,4>& bfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs,
			  TinyVector<double,4> &dbfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs,
			  TinyVector<complex<double>,4> &dbfuncs) const;
  inline void operator() (int i, TinyVector<double,4>& bfuncs,
			  TinyVector<double,4> &dbfuncs,
			  TinyVector<double,4> &d2bfuncs) const;
  inline void operator() (int i, TinyVector<complex<double>,4>& bfuncs,
			  TinyVector<complex<double>,4> &dbfuncs,
			  TinyVector<complex<double>,4> &d2bfuncs) const;
};

template<typename GridType> void
NUBsplineBasis<GridType>::Init(GridType *gridPtr, bool periodic)
{
  Periodic = periodic;
  GridPtr = gridPtr;
  GridType &grid =  (*GridPtr);

  int N = grid.NumPoints;
  dxInv.resize(grid.NumPoints+2);
  xVals.resize(N+4);
  for (int i=0; i<N; i++)
    xVals(i+2) = grid(i);
  if (!Periodic) {
    xVals(0) = grid(0) - 2.0*(grid(1)-grid(0));
    xVals(1) = grid(0) - 1.0*(grid(1)-grid(0));
    xVals(N+2) = grid(N-1) + 1.0*(grid(N-1)-grid(N-2));
    xVals(N+3) = grid(N-1) + 2.0*(grid(N-1)-grid(N-2));
  }
  else {
    xVals(1)   = grid(0)   - (grid(N-1) - grid(N-2));
    xVals(0)   = grid(0)   - (grid(N-1) - grid(N-3));
    xVals(N+2) = grid(N-1) + (grid(1)   - grid(0));
    xVals(N+3) = grid(N-1) + (grid(2)   - grid(0));
  }
  for (int i=0; i<N+2; i++) 
    for (int j=0; j<3; j++) 
      dxInv(i)[j] = 1.0/(xVals(i+j+1)-xVals(i));
}


template<typename GridType> int
NUBsplineBasis<GridType>::operator()(double x, TinyVector<double,4> &bfuncs) const
{
  GridType &grid =  (*GridPtr);
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
  return i;
}


template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<double,4> &bfuncs) const
{
  int i2 = i+2;
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  double x = grid(i);

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

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i,
				     TinyVector<complex<double>,4> &bfuncs) const
{
  TinyVector<double,4> funcs;
  (*this)(i, funcs);
  bfuncs[0] = complex<double>(funcs[0], funcs[0]);
  bfuncs[1] = complex<double>(funcs[1], funcs[1]);
  bfuncs[2] = complex<double>(funcs[2], funcs[2]);
  bfuncs[3] = complex<double>(funcs[3], funcs[3]);
}



template<typename GridType> int
NUBsplineBasis<GridType>::operator()(double x, TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs) const
{
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  int i = grid.ReverseMap (x);
  int i2 = i+2;

  b1[0]      = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]      = (x-xVals(i2))    * dxInv(i+2)[0];
  
  b2[0]      = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]      = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
		(xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]      = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  
  bfuncs[0]  = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1]  = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
		(xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2]  = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
		(xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3]  = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);

  return i;


  // The following is the stupid way of doing things:
  //   double  db1[2], db2[3];
  //   db1[0]     = -dxInv(i+2)[0];
  //   db1[1]     =  dxInv(i+2)[0];
  //   db2[0]     = ((xVals(i2+1)-x)  * dxInv(i+1)[1] * db1[0] - dxInv(i+1)[1] * b1[0]);
  //   db2[1]     = (((x-xVals(i2-1)) * dxInv(i+1)[1] * db1[0] + dxInv(i+1)[1] * b1[0]) +
  // 		((xVals(i2+2)-x) * dxInv(i+2)[1] * db1[1] - dxInv(i+2)[1] * b1[1]));
  //   db2[2]     = ((x-xVals(i2))    * dxInv(i+2)[1] * db1[1] + dxInv(i+2)[1] * b1[1]);
  //   dbfuncs[0] = ((xVals(i2+1)-x) * dxInv(i  )[2] * db2[0] - dxInv(i  )[2] * b2[0]);
  //   dbfuncs[1] = ((x-xVals(i2-2)) * dxInv(i  )[2] * db2[0] + dxInv(i  )[2] * b2[0] +
  // 		(xVals(i2+2)-x) * dxInv(i+1)[2] * db2[1] - dxInv(i+1)[2] * b2[1]);
  //   dbfuncs[2] = ((x-xVals(i2-1)) * dxInv(i+1)[2] * db2[1] + dxInv(i+1)[2] * b2[1] +
  // 		(xVals(i2+3)-x) * dxInv(i+2)[2] * db2[2] - dxInv(i+2)[2] * b2[2]);
  //   dbfuncs[3] = ((x-xVals(i2))   * dxInv(i+2)[2] * db2[2] + dxInv(i+2)[2] * b2[2]);
}



template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs) const
{
  int i2 = i+2;
  GridType &grid =  (*GridPtr);
  double x = grid(i);
  double b1[2], b2[3];

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

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);
}

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<complex<double>,4> &bfuncs,
				     TinyVector<complex<double>,4> &dbfuncs) const
{
  TinyVector<double,4> funcs, dfuncs;
  (*this)(i, funcs, dfuncs);
  bfuncs[0]  = complex<double>(funcs[0],  funcs[0]);
  bfuncs[1]  = complex<double>(funcs[1],  funcs[1]);
  bfuncs[2]  = complex<double>(funcs[2],  funcs[2]);
  bfuncs[3]  = complex<double>(funcs[3],  funcs[3]);
  dbfuncs[0] = complex<double>(dfuncs[0], dfuncs[0]);
  dbfuncs[1] = complex<double>(dfuncs[1], dfuncs[1]);
  dbfuncs[2] = complex<double>(dfuncs[2], dfuncs[2]);
  dbfuncs[3] = complex<double>(dfuncs[3], dfuncs[3]);
}

template<typename GridType> int
NUBsplineBasis<GridType>::operator()(double x, 
			    TinyVector<double,4> &bfuncs,
			    TinyVector<double,4> &dbfuncs,
			    TinyVector<double,4> &d2bfuncs) const
{
  GridType &grid =  (*GridPtr);
  double b1[2], b2[3];
  int i = min(grid.ReverseMap (x), grid.NumPoints-1);
  int i2 = i+2;

  b1[0]      = (xVals(i2+1)-x)  * dxInv(i+2)[0];
  b1[1]      = (x-xVals(i2))    * dxInv(i+2)[0];
  
  b2[0]      = (xVals(i2+1)-x)  * dxInv(i+1)[1] * b1[0];
  b2[1]      = ((x-xVals(i2-1)) * dxInv(i+1)[1] * b1[0]+
		(xVals(i2+2)-x) * dxInv(i+2)[1] * b1[1]);
  b2[2]      = (x-xVals(i2))    * dxInv(i+2)[1] * b1[1];
  
  bfuncs[0]  = (xVals(i2+1)-x)  * dxInv(i  )[2] * b2[0];
  bfuncs[1]  = ((x-xVals(i2-2)) * dxInv(i  )[2] * b2[0] +
		(xVals(i2+2)-x) * dxInv(i+1)[2] * b2[1]);
  bfuncs[2]  = ((x-xVals(i2-1)) * dxInv(i+1)[2] * b2[1] +
		(xVals(i2+3)-x) * dxInv(i+2)[2] * b2[2]);
  bfuncs[3]  = (x-xVals(i2))    * dxInv(i+2)[2] * b2[2];

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);

  d2bfuncs[0] = 6.0 * (+dxInv(i+0)[2]*dxInv(i+1)[1]*b1[0]);
  d2bfuncs[1] = 6.0 * (-dxInv(i+1)[1]*(dxInv(i+0)[2]+dxInv(i+1)[2])*b1[0] +
		       dxInv(i+1)[2]*dxInv(i+2)[1]*b1[1]);
  d2bfuncs[2] = 6.0 * (dxInv(i+1)[2]*dxInv(i+1)[1]*b1[0] -
		       dxInv(i+2)[1]*(dxInv(i+1)[2] + dxInv(i+2)[2])*b1[1]);
  d2bfuncs[3] = 6.0 * (+dxInv(i+2)[2]*dxInv(i+2)[1]*b1[1]);

  return i;
}

template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<complex<double>,4> &bfuncs,
				     TinyVector<complex<double>,4> &dbfuncs,
				     TinyVector<complex<double>,4> &d2bfuncs) const
{
  TinyVector<double,4> funcs, dfuncs, d2funcs;
  (*this)(i, funcs, dfuncs, d2funcs);
  bfuncs[0]   = complex<double>(  funcs[0],   funcs[0]);
  bfuncs[1]   = complex<double>(  funcs[1],   funcs[1]);
  bfuncs[2]   = complex<double>(  funcs[2],   funcs[2]);
  bfuncs[3]   = complex<double>(  funcs[3],   funcs[3]);
  dbfuncs[0]  = complex<double>( dfuncs[0],  dfuncs[0]);
  dbfuncs[1]  = complex<double>( dfuncs[1],  dfuncs[1]);
  dbfuncs[2]  = complex<double>( dfuncs[2],  dfuncs[2]);
  dbfuncs[3]  = complex<double>( dfuncs[3],  dfuncs[3]);
  d2bfuncs[0] = complex<double>(d2funcs[0], d2funcs[0]);
  d2bfuncs[1] = complex<double>(d2funcs[1], d2funcs[1]);
  d2bfuncs[2] = complex<double>(d2funcs[2], d2funcs[2]);
  d2bfuncs[3] = complex<double>(d2funcs[3], d2funcs[3]);
}



template<typename GridType> void
NUBsplineBasis<GridType>::operator()(int i, TinyVector<double,4> &bfuncs,
				     TinyVector<double,4> &dbfuncs,
				     TinyVector<double,4> &d2bfuncs) const
{
  int i2 = i+2;
  GridType &grid =  (*GridPtr);
  double x = grid(i);
  double b1[2], b2[3];

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

  dbfuncs[0] = -3.0 * (dxInv(i  )[2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv(i  )[2] * b2[0] - dxInv(i+1)[2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv(i+1)[2] * b2[1] - dxInv(i+2)[2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv(i+2)[2] * b2[2]);

  d2bfuncs[0] = 6.0 * (+dxInv(i+0)[2]*dxInv(i+1)[1]*b1[0]);
  d2bfuncs[1] = 6.0 * (-dxInv(i+1)[1]*(dxInv(i+0)[2]+dxInv(i+1)[2])*b1[0] +
		       dxInv(i+1)[2]*dxInv(i+2)[1]*b1[1]);
  d2bfuncs[2] = 6.0 * (dxInv(i+1)[2]*dxInv(i+1)[1]*b1[0] -
		       dxInv(i+2)[1]*(dxInv(i+1)[2] + dxInv(i+2)[2])*b1[1]);
  d2bfuncs[3] = 6.0 * (+dxInv(i+2)[2]*dxInv(i+2)[1]*b1[1]);
}



template<typename GridType, typename T> inline void
SolveDerivInterp1D (NUBsplineBasis<GridType> &basis,
		    Array<T,1> data, Array<T,1> p,
		    TinyVector<T,4> abcdInitial,
		    TinyVector<T,4> abcdFinal)
{
  assert (p.size() == (data.size()+2));

  // Banded matrix storage.  The first three elements in the
  // tinyvector store the tridiagonal coefficients.  The last element
  // stores the RHS data.
  Array<TinyVector<T,4>,1> bands(p.size());
  int M = data.size();

  // Fill up bands
  bands(0)          = abcdInitial;
  bands(p.size()-1) = abcdFinal;
  for (int i=0; i<M; i++) {
    basis (i, bands(i+1));
    bands(i+1)[3] = data(i);
  }
    
  // Now solve:
  // First and last rows are different
  bands(0)[1] /= bands(0)[0];
  bands(0)[2] /= bands(0)[0];
  bands(0)[3] /= bands(0)[0];
  bands(0)[0] = 1.0;
  bands(1)[1] -= bands(1)[0]*bands(0)[1];
  bands(1)[2] -= bands(1)[0]*bands(0)[2];
  bands(1)[3] -= bands(1)[0]*bands(0)[3];
  bands(0)[0] = 0.0;
  bands(1)[2] /= bands(1)[1];
  bands(1)[3] /= bands(1)[1];
  bands(1)[1] = 1.0;
  
  // Now do rows 2 through M+1
  for (int row=2; row < (p.size()-1); row++) {
    bands(row)[1] -= bands(row)[0]*bands(row-1)[2];
    bands(row)[3] -= bands(row)[0]*bands(row-1)[3];
    bands(row)[2] /= bands(row)[1];
    bands(row)[3] /= bands(row)[1];
    bands(row)[0] = 0.0;
    bands(row)[1] = 1.0;
  }

  // Do last row
  bands(M+1)[1] -= bands(M+1)[0]*bands(M-1)[2];
  bands(M+1)[3] -= bands(M+1)[0]*bands(M-1)[3];
  bands(M+1)[2] -= bands(M+1)[1]*bands(M)[2];
  bands(M+1)[3] -= bands(M+1)[1]*bands(M)[3];
  bands(M+1)[3] /= bands(M+1)[2];
  bands(M+1)[2] = 1.0;

  p(M+1) = bands(M+1)[3];
  // Now back substitute up
  for (int row=M; row>0; row--)
    p(row) = bands(row)[3] - bands(row)[2]*p(row+1);
  
  // Finish with first row
  p(0) = bands(0)[3] - bands(0)[1]*p(1) - bands(0)[2]*p(2);
}

template<typename T, int N> inline TinyVector<T,4> 
real(TinyVector<complex<T>,N> v)
{
  TinyVector<T,N> rv;
  for (int i=0; i<N; i++)
    rv[i] = real(v[i]);
  return rv;
}

template<typename T, int N> inline TinyVector<T,4> 
imag(TinyVector<complex<T>,N> v)
{
  TinyVector<T,N> iv;
  for (int i=0; i<N; i++)
    iv[i] = imag(v[i]);
  return iv;
}


template<typename GridType, typename T> inline void
SolveDerivInterp1D (NUBsplineBasis<GridType> &basis,
		    Array<complex<T>,1> data, Array<complex<T>,1> p,
		    TinyVector<complex<T>,4> abcdInitial,
		    TinyVector<complex<T>,4> abcdFinal)
{
  Array<T,1> rdata(data.size()), rp(p.size()),
    idata(data.size()), ip(p.size());
  TinyVector<T,4> rStartBC, iStartBC, rEndBC, iEndBC;
  rStartBC = real(abcdInitial);
  iStartBC = imag(abcdInitial);
  rEndBC   = real(abcdFinal);
  iEndBC   = imag(abcdFinal);
  
  for (int i=0; i<data.size(); i++){
    rdata(i) = real(data(i));
    idata(i) = imag(data(i));
  }
  SolveDerivInterp1D (basis, rdata, rp, rStartBC, rEndBC);
  SolveDerivInterp1D (basis, idata, ip, iStartBC, iEndBC);
  for (int i=0; i<p.size(); i++)
    p(i) = complex<T> (rp(i), ip(i));
}


template<typename GridType, typename T> inline void
SolvePeriodicInterp1D (NUBsplineBasis<GridType> &basis,
		       Array<T,1> data, Array<T,1> p)
{
  assert (p.size() == (data.size()+3));

  // Banded matrix storage.  The first three elements in the
  // tinyvector store the tridiagonal coefficients.  The last element
  // stores the RHS data.
  Array<TinyVector<double,4>,1> bands(data.size());
  Array<double,1> lastCol(data.size());
  int M = data.size();

  // Fill up bands
  for (int i=0; i<M; i++) {
    basis (i, bands(i));
    bands(i)[3] = data(i);
  }
    
  // Now solve:
  // First and last rows are different
  bands(0)[2] /= bands(0)[1];
  bands(0)[0] /= bands(0)[1];
  bands(0)[3] /= bands(0)[1];
  bands(0)[1]  = 1.0;
  bands(M-1)[1] -= bands(M-1)[2]*bands(0)[0];
  bands(M-1)[3] -= bands(M-1)[2]*bands(0)[3];
  bands(M-1)[2]  = -bands(M-1)[2]*bands(0)[2];
  lastCol(0) = bands(0)[0];
  
  for (int row=1; row < (M-1); row++) {
    bands(row)[1] -= bands(row)[0] * bands(row-1)[2];
    bands(row)[3] -= bands(row)[0] * bands(row-1)[3];
    lastCol(row)   = -bands(row)[0] * lastCol(row-1);
    bands(row)[0] = 0.0;
    bands(row)[2] /= bands(row)[1];
    bands(row)[3] /= bands(row)[1];
    lastCol(row)  /= bands(row)[1];
    bands(row)[1]  = 1.0;
    if (row < (M-2)) {
      bands(M-1)[3] -= bands(M-1)[2]*bands(row)[3];
      bands(M-1)[1] -= bands(M-1)[2]*lastCol(row);
      bands(M-1)[2] = -bands(M-1)[2]*bands(row)[2];
    }
  }

  // Now do last row
  // The [2] element and [0] element are now on top of each other 
  bands(M-1)[0] += bands(M-1)[2];
  bands(M-1)[1] -= bands(M-1)[0] * (bands(M-2)[2]+lastCol(M-2));
  bands(M-1)[3] -= bands(M-1)[0] *  bands(M-2)[3];
  bands(M-1)[3] /= bands(M-1)[1];
  p(M) = bands(M-1)[3];
  for (int row=M-2; row>=0; row--) 
    p(row+1) = bands(row)[3] - bands(row)[2]*p(row+2) - lastCol(row)*p(M);
  
  p(0) = p(M);
  p(M+1) = p(1);
  p(M+2) = p(2);
}


template<typename GridType, typename T> inline void
SolvePeriodicInterp1D (NUBsplineBasis<GridType> &basis,
		       Array<complex<T>,1> data, Array<complex<T>,1> p)
{
  Array<T,1> rdata(data.size()), rp(p.size()),
    idata(data.size()), ip(p.size());
  
  for (int i=0; i<data.size(); i++){
    rdata(i) = real(data(i));
    idata(i) = imag(data(i));
  }
  SolvePeriodicInterp1D (basis, rdata, rp);
  SolvePeriodicInterp1D (basis, idata, ip);
  for (int i=0; i<p.size(); i++)
    p(i) = complex<T> (rp(i), ip(i));
}



#endif
