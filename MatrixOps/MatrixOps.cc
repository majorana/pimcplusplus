#include "MatrixOps.h"

extern "C" void dgesvd_(char *JOBU, char* JOBVT, int *M, int *N,
			double *A, int *LDA, double *S, double *U,
			int *LDU, double *VT, int *LDVT, double *work,
			int *LWORK, int *INFO);


void SVdecomp (Array<double,2> &A,
	       Array<double,2> &U, Array<double,1> &S,
	       Array<double,2> &V)
{
  int M = A.rows();
  int N = A.cols();
  cerr << " M = " << M << " N = " << N << endl;
  Array<double,2> Atrans(M,N);
  // U will be Utrans after lapack call
  U.resize(M,M);
  V.resize(N,N);
  
  S.resize(min(N,M));
  Atrans = A;

  // Transpose U for FORTRAN ordering
  Transpose(Atrans);
  char JOBU = 'A';   // return all columns of U
  char JOBVT = 'A'; // return all columsn of V
  int LDA = M;
  int LDU = M;
  int LDVT = N;
  int LWORK = 10 * max(3*min(M,N)+max(M,N),5*min(M,N));
  Array<double,1> WORK(LWORK);
  int INFO;

  dgesvd_(&JOBU, &JOBVT, &M, &N, Atrans.data(), &LDA,
	  S.data(), U.data(), &LDU, V.data(), &LDVT,
	  WORK.data(), &LWORK, &INFO);
  assert (INFO == 0);
  // Transpose U to get back to C ordering
  // V was really Vtrans so we don't need to transpose
  Transpose(U);
  Array<double,2> tmp(M,min(N,M));
  for (int i=0; i<M; i++)
    for (int j=0; j<min(N,M);j++)
      tmp(i,j) = U(i,j);
  U.resize(M, min(N,M));
  U = tmp;
  cerr << "S = " << S << endl;
}



	       

      // Adapted from Numerical Recipes in C
void LUdecomp (Array<double,2> &A, Array<int,1> &perm, 
	       double &sign)
{
  int i, imax, j, k;
  int n = A.rows();
  double big, dum, sum, temp;
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
	sum -= A(i,k)*A(k,j);
      A(i,j) = sum;
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
      dum=1.0/A(j,j);
      for (i=j+1; i<n; i++)
	A(i,j) *= dum;
    }
  }  
}


void LUsolve (Array<double,2> &LU, Array<int,1> &perm,
	      Array<double,1> &b)
{
  int i, ii=-1,ip,j;
  double sum;
  int n = LU.rows();
  
  for (i=0; i<n; i++) {
    ip = perm(i);
    sum = b(ip);
    b(ip) = b(i);
    if (ii>=0)
      for (j=ii; j<i; j++)
	sum -= LU(i,j)*b(j);
    else if (sum) 
      ii = i;
    b(i) = sum;
  }
  for (i=n-1; i>=0; i--) {
    sum = b(i);
    for (j=i+1; j<n; j++)
      sum -= LU(i,j)*b(j);
    b(i) = sum/LU(i,i);
  }
}


Array<double,2> Inverse (Array<double,2> &A)
{
  Array<double,2> LU(A.rows(), A.cols()), Ainv(A.rows(), A.cols());
  Array<double,1> col(A.rows());
  Array<int,1> perm;
  double sign;

  LU = A;
  LUdecomp (LU, perm, sign);
  for (int j=0; j<A.rows(); j++) {
    for (int i=0; i<A.rows(); i++)
      col(i) = 0.0;
    col(j) = 1.0;
    LUsolve (LU, perm, col);
    for (int i=0; i<A.rows(); i++)
      Ainv(i,j) = col(i);
  }
  return (Ainv);
}

  
