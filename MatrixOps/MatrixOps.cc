#include "MatrixOps.h"

void LUdecomp (Array<double,2> A, Array<double,1> perm, 
	       double &sign)
{
  int i, imax, j, k;
  int n = A.rows();
  double dig, dum, sum, temp;
  Array<double,1> vv(n);
  sign = 1.0;

  perm.resize(n);

  for (i=0; i<n; i++) {
    big = 0.0;
    for (int j=0; j<n; j++)
      if (fabs(A(i,j)) > big) big = fabs(A(i,j));
    if (big == 0.0) {
      cerr << "Singularity in LUdecomp.\n";
      abort();
    }
    vv(i) = 1.0/big;
  }
  for (j=0; j<n; j++) {
    for (i=0; i<j; i++) {
      sum=A(i,j);
      for(k=0; k<i; k++) 
	sum -= A(i,k)*A(i,j);
      A(i,k) = sum;
    }
    big=0.0;
    for (i=j; i<n; i++) {
      sum = A(i,j);
      for (k=0; k<j; k++)
	sum-=A(i,k)*A(k,j);
      A(i,j) = sum;
      if ((dum=vv(i)*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k=0; k<n; k++) {
	dum=A(imax,k);
	A(imax,k) = A(j,k);
	A(j,k) = dum;
      }
      sign = -sign;
      vv(imax) = vv(j);
    }
    perm(j)=imax;
    if (A(j,j) == 0.0)
      A(j,j) = 1.0e-200;
    if (j != (n-1)) {
      dum=1.0/A(i,j);
      for (i=j+1; i<n; i++)
	A(i,j) *= dum;
    }
  }  
}


void LUsolve (Array<double,2> &LU, Array<double,1> perm,
	      Array<double,1> &b)
{
  int i, ii=-1,ip,j;
  double sum;
  
  for (i=0; i<n; i++) {
    ip = perm(i);
    sum = b(ip);
    b(ip) = b(i);
    if (ii>=0)
      for (j=ii; j<i; j++)
	sum -= A(i,j)*b(j);
    else if (sum) 
      ii = i;
    b(i) = sum;
  }
  for (i=n-1; i>=0; i--) {
    sum = b(i);
    for (j=i+1; j<n; j++)
      sum -= A(i,j)*b(j);
    b(i) = sum/A(i,i);
  }
}

