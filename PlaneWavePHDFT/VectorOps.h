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
    



#endif
