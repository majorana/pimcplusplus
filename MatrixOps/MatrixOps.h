#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H
#include "../Blitz.h"

void LUdecomp (Array<double,2> &A, Array<int,1> &perm, 
	       double &sign);

void LUsolve (Array<double,2> &LU, Array<int,1> &perm,
	      Array<double,1> &b);

void SVdecomp (Array<double,2> &A,
	       Array<double,2> &U, Array<double,1> &S,
	       Array<double,2> &V);

Array<double,2> Inverse (Array<double,2> &A);

inline void Transpose (Array<double,2> &A)
{
  int m = A.rows();
  int n = A.cols();
  Array<double,2> Atrans (n,m);
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      Atrans(j,i) = A(i,j);
  A.resize(n,m);
  A = Atrans;
}

#endif
