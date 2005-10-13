#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

// extern "C"{
// #ifdef USE_MKL
//   #include <mkl_cblas.h>
// #else  
//   #include <cblas.h> 
// #endif
// }
#include "../Blitz.h"
#include "../config.h"

typedef TinyVector<int,3> Int3;
typedef Array<complex<double>,1> zVec;
typedef Array<cVec3,1> zVecVec;
typedef Array<cMat3,1> zMatVec;

#define F77_DZNRM2 F77_FUNC(dznrm2,DZNRM2)
#define F77_ZDSCAL F77_FUNC(zdscal,ZDSCAL)
#define F77_ZDOTC  F77_FUNC(zdotc,ZDOTC)
#define F77_ZGEMV  F77_FUNC(zgemv,ZGEMV)
#define F77_ZAXPY  F77_FUNC(zaxpy,ZAXPY)

extern "C" double F77_DZNRM2(const int *N, const void *X, const int *INC);
extern "C" void   F77_ZDSCAL(const int *N, double *ALPHA, const void *X, 
			     const int *INC);
extern "C" complex<double> F77_ZDOTC (const int *N, 
				      const void *X, const int *INCX, 
				      const void *Y, const int *INCY);
extern "C" void   F77_ZGEMV (char *TRANS, const int *M, const int *N, 
			     complex<double> *alpha, const void *A, 
			     const int *LDA, const void *X, 
			     const int *INCX, complex<double> *beta, 
			     const void *Y, const int *INCY);

extern "C" void   F77_ZAXPY (const int *N, complex<double> *ALPHA,
			     void *X, int *INCX, void *Y, int *INCY);

inline void Normalize (zVec &c)
{
//   double norm = 0.0;
//   for (int i=0; i<c.size(); i++)
//     norm += c(i).real()*c(i).real() + c(i).imag()*c(i).imag();
//   cerr << "norm1 = " << norm << endl;
//   norm = 1.0/sqrt(norm);
//   for (int i=0; i<c.size(); i++)
//     c(i) *= norm;

//   double norm = cblas_dznrm2(c.size(), c.data(), 1);
//   cblas_zdscal(c.size(), 1.0/norm, c.data(), 1);

  const int inc=1;
  int n = c.size();
  double norm = F77_DZNRM2(&n, c.data(), &inc);
  norm = 1.0/norm;
  F77_ZDSCAL(&n, &norm, c.data(), &inc);
}

inline double norm (const zVec &c)
{
  double norm;
  int n = c.size();
  const int inc=1;
  return F77_DZNRM2(&n, c.data(), &inc);
    //  return cblas_dznrm2(c.size(), c.data(), 1);
}

inline complex<double> conjdot(zVec &cA, zVec &cB)
{
 //  complex<double> z(0.0, 0.0);
//   for (int i=0; i<cA.size(); i++)
//     z += conj(cA(i))*cB(i);
//   return z;
  const int n = cA.size();
  const int inc = 1;
  complex<double> z;
  return F77_ZDOTC (&n, cA.data(), &inc, cB.data(), &inc);
//   cblas_zdotc_sub (cA.size(), cA.data(), 1, cB.data(), 1, &z);
//  return z;
}

inline double realconjdot(zVec &cA, zVec &cB)
{
//   double re = 0.0;
//   for (int i=0; i<cA.size(); i++)
//     re += real(conj(cA(i))*cB(i));
//   return re;

  // complex<double> z;
//   cblas_zdotc_sub (cA.size(), cA.data(), 1, cB.data(), 1, &z);
//   return z.real();
  return conjdot(cA,cB).real();
}

// inline Array<complex<double>,3>& operator*=
// (Array<complex<double>,3> &a, Array<complex<double>,3> &b)
// {
//   int N = a.size();
//   complex<double> *aPtr, *bPtr;
//   aPtr = a.data();
//   bPtr = b.data();
//   for (int i=0; i<N; i++) {
//     *aPtr *= *bPtr;
//     aPtr++;
//     bPtr++;
//   }
// }
    

inline void 
Orthogonalize (const Array<complex<double>,2> &A, zVec &x)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  complex<double> zero(0.0, 0.0);
  complex<double> one (1.0, 0.0);
  complex<double> minusone (-1.0, 0.0);
  Array<complex<double>,1> S(m);
  
  // Calculate overlaps
  // Calling with column major and ConjTrans is equivalent to
  // conjugate of untransposed row major
  char trans='C';
  const int inc=1;

  F77_ZGEMV(&trans, &n, &m, &one, A.data(), &n, x.data(), &inc,
	    &zero, S.data(), &inc);

//   cblas_zgemv(CblasColMajor, CblasConjTrans, n, m, &one,
// 	      A.data(), n, x.data(), 1, &zero, S.data(), 1);

//   for (int i=0; i<m; i++) {
//     fprintf (stderr, "S[%d] = %18.14f + %18.14fi\n",
// 	     real(S[i]), imag(S[i]));
    
  // Now, subtract off components * overlaps
  trans='T';
  F77_ZGEMV(&trans, &m, &n, &minusone, A.data(), &n,
	    S.data(), &inc, &one, x.data(), &inc);
//   cblas_zgemv(CblasRowMajor, CblasTrans, m, n, &minusone,
//  	      A.data(), n, S.data(), 1, &one, x.data(), 1);

}

inline double mag (complex<double> x)
{
  return (x.real()*x.real() + x.imag()*x.imag());
}

inline void 
Orthogonalize2 (Array<complex<double>,2> &A, zVec &x, int lastBand)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;
  Array<complex<double>,1> S(m);

  for (int row=0; row<=lastBand; row++) {
    Ar.reference (A(row,Range::all()));
    S(row) = conjdot (Ar, x);
  }
  for (int row=0; row<=lastBand; row++) 
    x -= S(row) * A(row,Range::all());
//   for (int row=0; row<=lastBand; row++) {
//     Ar.reference (A(row,Range::all()));
//     S(row) = conjdot (Ar, x);
//     if (mag(S(row)) > 1.0e-14) {
//       cerr << "row = " << row << " lastband = " << lastBand << endl;
//       cerr << "Error in Orthogonalize2!, s = " << S(row) << endl;
//       double norm = realconjdot (Ar, Ar);
//       cerr << "norm = " << norm << endl;
//     }
//   }
}

inline void
OrthogExcluding(const Array<complex<double>,2> &A, zVec &x,
		int excluding)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;
  complex<double> S[m];

  for (int row=0; row<m; row++) {
    Ar.reference (A(row,Range::all()));
    S[row] = conjdot (Ar, x);
  }
  for (int row=0; row<m; row++) 
    if (row != excluding)
      x -= S[row] * A(row,Range::all());
}


inline void
OrthogLower(const Array<complex<double>,2> &A, zVec &x,
	    int currBand)
{
  int m = currBand;
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;
  complex<double> S;

  for (int row=0; row<m; row++) {
    Ar.reference (A(row,Range::all()));
    S = conjdot (Ar, x);
    x -= S * A(row,Range::all());
  }
}


inline void
GramSchmidt (Array<complex<double>,2> &A)
{
  zVec a, b;
  for (int i=0; i<A.rows(); i++) {
    a.reference (A(i,Range::all()));
    Normalize(a);
    for (int j=i+1; j<A.rows(); j++) {
      b.reference(A(j,Range::all()));
      b = b - (conjdot(a,b)*a);
      Normalize(b);
    }
  }
}

inline void
Orthogonalize (Array<complex<double>,2> &A)
{
  zVec x, y;
  Array<complex<double>,1> S(A.rows());
  for (int iter=0; iter < 40; iter++) {
    for (int i=0; i<A.rows();i++) {
      x.reference (A(i, Range::all()));
      for (int j=i+1; j<A.rows(); j++) {
	y.reference (A(j,Range::all()));
	S(j) = conjdot (y, x);
      }
      for (int j=i+1; j<A.rows(); j++) {
      y.reference (A(j,Range::all()));
      x -= S(j) * y;
      }
      Normalize (x);
    }
  }
}


inline void CheckOrthog (const Array<complex<double>,2> &A,
			 zVec &x)
{
  zVec Ai;
  double normInv = 1.0/norm(x);
  for (int i=0; i<A.rows(); i++) {
    Ai.reference (A(i,Range::all()));
    if (normInv*mag(conjdot(Ai, x)) > 1.0e-13) {
      cerr << "CheckOrthog failed for i=" << i << ".\n";
      exit(1);
    }
  }
}

inline void
zaxpy (complex<double> alpha, const Array<complex<double>,1> &x,
       const complex<double> y, Array<complex<double>,1> &axpy)
{

}


// inline Array<complex<double>,3>&
// operator*= (Array<complex<double>,3> &A, const Array<complex<double>,3> &B)
// {
  

// }


#endif
