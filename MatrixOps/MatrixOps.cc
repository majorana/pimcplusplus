#include "MatrixOps.h"
#include <blitz/mstruct.h>

extern "C" void dgesvd_(char *JOBU, char* JOBVT, int *M, int *N,
			double *A, int *LDA, double *S, double *U,
			int *LDU, double *VT, int *LDVT, double *work,
			int *LWORK, int *INFO);

extern "C" void dgetrf_(int *m, int *n, double A[], int *lda, int ipiv[],
			int *info);

extern "C" void dgemm_(char *transA, char *transB, int *m, int *n, int *k,
		       double *alpha, const double *A, int *lda, const double *B, int *ldb,
		       double *beta,  double *C, int *ldc);

const Array<double,2> operator*(const Array<double,2> &A,
				const Array<double,2> &B)
{
  int m = A.rows();
  int n = B.cols();
  int k = A.cols();
  assert (B.rows() == k);
  // We use "transpose" operation because we have C ordering, which fortran
  // thinks is transposed.
  char transA = 'T';
  char transB = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  GeneralArrayStorage<2> colMajor;
  colMajor.ordering() = firstDim, secondDim;
  Array<double,2> C(m,n,colMajor);
  dgemm_(&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
	 B.data(), &n,
	 &beta, C.data(), &m);
  return C;
}


void MatMult (const Array<double,2> &A, const Array<double,2> &B,
	      Array<double,2> &C)
{
  int m = A.rows();
  int n = B.cols();
  int k = A.cols();
  assert (B.rows() == k);
  // We use "transpose" operation because we have C ordering, which fortran
  // thinks is transposed.
  char transA = 'T';
  char transB = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  GeneralArrayStorage<2> colMajor;
  colMajor.ordering() = firstDim, secondDim;
  dgemm_(&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
	 B.data(), &n,
	 &beta, C.data(), &m);
}


double Determinant (const Array<double,2> &A)
{
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take a determinant of a non-square
		    // matrix
  if (A.rows() == 2) 
    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
  else {
    Array<double,2> LU(m,m);
    Array<int,1> ipiv(m);
    int info;
    LU = A;
    // Do LU factorization
    dgetrf_(&m, &n, LU.data(), &m, ipiv.data(), &info);
    double det = 1.0;
    int numPerm = 0;
    for (int i=0; i<m; i++) {
      det *= LU(i,i);
      numPerm += (ipiv(i) != (i+1));
    }
    if (numPerm & 1)
      det *= -1.0;
    
    return det;
  }
}

// Replaces A with its inverse by gauss-jordan elimination with full pivoting
// Adapted from Numerical Recipes in C
void GJInverse (Array<double,2> &A)
{
  assert (A.cols() == A.rows());
  int n = A.rows();
  int colIndex[n], rowIndex[n], ipiv[n];
  double big, dum, pivInv, temp;
  int icol, irow;
  
  for (int j=0; j<n; j++)
    ipiv[j] = -1;

  for (int i=0; i<n; i++) {
    big = 0.0;
    for (int j=0; j<n; j++) 
      if (ipiv[j] != 0)
	for (int k=0; k<n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs(A(j,k)) >= big) {
	      big = fabs(A(j,k));
	      irow = j; 
	      icol = k;
	    }
	  }
	  else if (ipiv[k] > 0) {
	    cerr << "GJInverse: Singular matrix!\n";
	    abort();
	  }
	}
    ++(ipiv[icol]); 
    
    if (irow != icol) 
      for (int l=0; l<n; l++) 
	swap (A(irow,l), A(icol,l));
    
    rowIndex[i] = irow;
    colIndex[i] = icol;
    if (A(icol,icol) == 0.0) { 
      cerr << "GJInverse: Singular matrix!\n";
      abort();
    }
    pivInv = 1.0/A(icol,icol);
    A(icol,icol) = 1.0;
    for (int l=0; l<n; l++)
      A(icol,l) *= pivInv;
    for (int ll=0; ll<n; ll++)
      if (ll != icol) {
	double dum = A(ll,icol);
	A(ll,icol) = 0.0;
	for (int l=0; l<n; l++)
	  A(ll,l) -= A(icol,l)*dum;
      }
  }
  // Now unscramble the permutations
  for (int l=n-1; l>=0; l--) {
    if (rowIndex[l] != colIndex[l])
      for (int k=0; k<n ; k++)
	swap (A(k,rowIndex[l]),A(k, colIndex[l]));
  }
}


inline void SwapRow (Array<double,2> A, int row1, int row2)
{
  int m = A.cols();
  for (int col=0; col<m; col++) {
    double temp = A(row1,col);
    A(row1,col) = A(row2,col);
    A(row2,col) = temp;
  }
}


// The cofactors of A are given by 
// cof(A) = det(A) transpose(A^{-1})
void Cofactors (const Array<double,2> &A, 
		Array<double,2> &cof,
		Array<double,2> &scratch)
{
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take cofactors of a non-square matrix
  int  ipiv[m];
  int info;
  // Copy and transpose for FORTRAN ordering
  for (int i=0; i<m; i++)
    for (int j=0; j<m; j++)
      scratch(i,j) = A(j,i);
  // Do LU decomposition
  dgetrf_ (&m, &n, scratch.data(), &m, ipiv, &info);
  // Now scratch contains LU matrix in fortran ordering with pivots in ipiv
  // Put identity matrix in cof
  cof = 0.0;
  for (int i=0; i<m; i++)
    cof(i,i) = 1.0;
  int numPerm = 0;
  // Now apply permutation matrix to cof
  for (int row=0; row<m; row++) {
    int ip = ipiv[row]-1;
    if (ip != row) { 
      SwapRow(cof, row, ip);
      numPerm++;
    }
  }
}


void SVdecomp (Array<double,2> &A,
	       Array<double,2> &U, Array<double,1> &S,
	       Array<double,2> &V)
{
  int M = A.rows();
  int N = A.cols();
  Array<double,2> Atrans(M,N);
  // U will be Utrans after lapack call
  U.resize(min(M,N),M);
  V.resize(N,min(M,N));
  
  S.resize(min(N,M));
  Atrans = A;

  // Transpose U for FORTRAN ordering
  Transpose(Atrans);
  char JOBU  = 'S'; // return min (M,N) columns of U
  char JOBVT = 'S'; // return min (M,N) columns of V
  int LDA = M;
  int LDU = M;
  int LDVT = min(M,N);
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

  
