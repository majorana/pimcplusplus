#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

extern "C"{  
#include <cblas.h> 
}
#include "../Blitz.h"

typedef TinyVector<int,3> Int3;
typedef Array<complex<double>,1> zVec;
typedef Array<cVec3,1> zVecVec;
typedef Array<cMat3,1> zMatVec;

inline void Normalize (zVec &c)
{
//   double norm = 0.0;
//   for (int i=0; i<c.size(); i++)
//     norm += c(i).real()*c(i).real() + c(i).imag()*c(i).imag();
//   cerr << "norm1 = " << norm << endl;
//   norm = 1.0/sqrt(norm);
//   for (int i=0; i<c.size(); i++)
//     c(i) *= norm;
  double norm = cblas_dznrm2(c.size(), c.data(), 1);
  cblas_zdscal(c.size(), 1.0/norm, c.data(), 1);
}

inline complex<double> conjdot(zVec &cA, zVec &cB)
{
 //  complex<double> z(0.0, 0.0);
//   for (int i=0; i<cA.size(); i++)
//     z += conj(cA(i))*cB(i);
//   return z;
  complex<double> z;
  cblas_zdotc_sub (cA.size(), cA.data(), 1, cB.data(), 1, &z);
  return z;
}

inline double realconjdot(zVec &cA, zVec &cB)
{
//   double re = 0.0;
//   for (int i=0; i<cA.size(); i++)
//     re += real(conj(cA(i))*cB(i));
//   return re;
  complex<double> z;
  cblas_zdotc_sub (cA.size(), cA.data(), 1, cB.data(), 1, &z);
  return z.real();
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
  complex<double> S[m];
  
  // Calculate overlaps
  // Calling with column major and ConjTrans is equivalent to
  // conjugate of untransposed row major
  cblas_zgemv(CblasColMajor, CblasConjTrans, n, m, &one,
	      A.data(), n, x.data(), 1, &zero, S, 1);

//   for (int i=0; i<m; i++) {
//     fprintf (stderr, "S[%d] = %18.14f + %18.14fi\n",
// 	     real(S[i]), imag(S[i]));
    
  // Now, subtract off components * overlaps
  cblas_zgemv(CblasRowMajor, CblasTrans, m, n, &minusone,
 	      A.data(), n, S, 1, &one, x.data(), 1);

}

inline double mag (complex<double> x)
{
  return (x.real()*x.real() + x.imag()*x.imag());
}

inline void
Orthogonalize (Array<complex<double>,2> &A);

inline void 
Orthogonalize2 (Array<complex<double>,2> &A, zVec &x,
		int exclBand)
{
  int m = A.rows();
  int n = A.cols();
  assert (n == x.size());
  zVec Ar;
  complex<double> S[m];

  for (int row=0; row<A.rows(); row++) {
    Ar.reference (A(row,Range::all()));
    S[row] = conjdot (Ar, x);
  }
  for (int row=0; row<A.rows(); row++) 
    if (row != exclBand)
      x -= S[row] * A(row,Range::all());
  for (int row=0; row<A.rows(); row++) {
    Ar.reference (A(row,Range::all()));
    S[row] = conjdot (Ar, x);
    if ((row != exclBand) && (mag(S[row]) > 1.0e-14))
      cerr << "Error in Orthogonalize2!, s = " << S[row] << endl;
  }
}

inline void
Orthogonalize (Array<complex<double>,2> &A)
{
  zVec x, y;
  complex<double> S[A.rows()];
  for (int iter=0; iter < 40; iter++) {
    for (int i=0; i<A.rows();i++) {
      x.reference (A(i, Range::all()));
      for (int j=i+1; j<A.rows(); j++) {
	y.reference (A(j,Range::all()));
	S[j] = conjdot (y, x);
      }
      for (int j=i+1; j<A.rows(); j++) {
      y.reference (A(j,Range::all()));
      x -= S[j] * y;
      }
      Normalize (x);
    }
  }
}


inline void CheckOrthog (const Array<complex<double>,2> &A,
			 zVec &x)
{
  zVec Ai;
  for (int i=0; i<A.rows(); i++) {
    Ai.reference (A(i,Range::all()));
    if (mag(conjdot(Ai, x)) > 1.0e-14) {
      cerr << "CheckOrthog failed.\n";
      exit(1);
    }
  }
}


#endif
