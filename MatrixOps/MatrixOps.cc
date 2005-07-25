#include "MatrixOps.h"
#include <blitz/mstruct.h>
#include "../config.h"

//#ifdef NOUNDERSCORE 
//#define FORT(name) name
//#else
#define FORT(name) name ## _
//#endif 


#define F77_DGESVD F77_FUNC(dgesvd,DGESVD)
#define F77_DGETRF F77_FUNC(dgetrf,DGETRF)
#define F77_ZGETRF F77_FUNC(zgetrf,ZGETRF)
#define F77_DGETRI F77_FUNC(dgetri,DGETRI)
#define F77_ZGETRI F77_FUNC(zgetri,ZGETRI)
#define F77_DGEMM  F77_FUNC(dgemm,DGEMM)
#define F77_DSYEVR F77_FUNC(dsyevr,DSYEVR)
#define F77_ZHEEVR F77_FUNC(zheevr,ZHEEVR)

extern "C" void 
F77_DGESVD (char *JOBU, char* JOBVT, int *M, int *N,
	    double *A, int *LDA, double *S, double *U,
	    int *LDU, double *VT, int *LDVT, double *work,
	    int *LWORK, int *INFO);

extern "C" void 
F77_DGETRF(int *m, int *n, double A[], int *lda, int ipiv[], int *info);

extern "C" void 
F77_ZGETRF(int *m, int *n, complex<double> A[], 
	   int *lda, int ipiv[], int *info);

extern "C" void 
F77_DGETRI (int *N, double A[], int *lda, int ipiv[], double work[], 
	    int *lwork, int *info);

extern "C" void 
F77_ZGETRI (int *N, complex<double> A[], int *lda, int ipiv[], 
	    complex<double> work[], int *lwork, int *info);

extern "C" void 
F77_DGEMM (char *transA, char *transB, int *m, int *n, int *k,
	   double *alpha, const double *A, int *lda, const double *B, int *ldb,
	   double *beta,  double *C, int *ldc);

extern "C" void 
F77_DSYEVR (char *JobType, char *Range, char *UpperLower, 
	    int *N, double *Amat, int *LDA,
	    double *VL, double *VU,
	    int *IL, int *IU, 
	    double *AbsTolerance, int *M,
	    double *EigVals, 
	    double *EigVecs, int *LDEigVecs, int *ISuppZ,
	    double *Work, int *Lwork, 
	    int *IWorkSpace, int *LIwork,
	    int *Info);

extern "C" void 
F77_ZHEEVR (char *JobType, char *Range, char *UpperLower, 
	    int *N, complex<double> *Amat, int *LDA,
	    double *VL, double *VU,
	    int *IL, int *IU, 
	    double *AbsTolerance, int *M,
	    double *EigVals, 
	    complex<double> *EigVecs, int *LDEigVecs, 
	    int *ISuppZ,
	    complex<double> *Work, int *Lwork, 
	    double *Rwork, int *LRwork,
	    int *IWorkSpace, int *LIwork,
	    int *Info);


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
  F77_DGEMM (&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
	     B.data(), &n, &beta, C.data(), &m);
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
  F77_DGEMM (&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
	     B.data(), &n, &beta, C.data(), &m);
}


double Determinant (const Array<double,2> &A)
{
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take a determinant of a non-square
		    // matrix
  if (A.rows() == 1)
    return (A(0,0));
  if (A.rows() == 2) 
    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
  else {
    Array<double,2> LU(m,m);
    Array<int,1> ipiv(m);
    int info;
    LU = A;
    // Do LU factorization
    F77_DGETRF (&m, &n, LU.data(), &m, ipiv.data(), &info);
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

complex<double> 
Determinant (const Array<complex<double>,2> &A)
{
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take a determinant of a non-square
		    // matrix
  if (A.rows() == 1)
    return (A(0,0));
  if (A.rows() == 2) 
    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
  else {
    Array<complex<double>,2> LU(m,m);
    Array<int,1> ipiv(m);
    int info;
    LU = A;
    // Do LU factorization
    F77_ZGETRF (&m, &n, LU.data(), &m, ipiv.data(), &info);
    complex<double> det = 1.0;
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
  const int maxSize = 2000;
  assert (A.cols() == A.rows());
  assert (A.cols() <= maxSize);
  int n = A.rows();

  if (n == 2) { // Special case for 2x2
    double a=A(0,0); double b=A(0,1);
    double c=A(1,0); double d=A(1,1);
    double detInv = 1.0/(a*d-b*c);
    A(0,0) = d*detInv;
    A(0,1) = -b*detInv;
    A(1,0) = -c*detInv;
    A(1,1) = a*detInv;
    return;
  }

  int colIndex[maxSize], rowIndex[maxSize], ipiv[maxSize];
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
	    cerr << "A = " << A << endl;
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
      cerr << "A = " << A << endl;
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
  const int maxSize = 2000;
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take cofactors of a non-square matrix
  assert (A.cols() < maxSize);
  int  ipiv[maxSize];
  int info;
  // Copy and transpose for FORTRAN ordering
  for (int i=0; i<m; i++)
    for (int j=0; j<m; j++)
      scratch(i,j) = A(j,i);
  // Do LU decomposition
  F77_DGETRF (&m, &n, scratch.data(), &m, ipiv, &info);
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

  F77_DGESVD (&JOBU, &JOBVT, &M, &N, Atrans.data(), &LDA,
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

  
void SymmEigenPairs (const Array<scalar,2> &A, int NumPairs,
		     Array<scalar,1> &Vals,
		     Array<scalar,2> &Vectors)
{
  char JobType = 'V';    // Find eigenvectors and eignevalues
  char Range   = 'I';    // Find eigenpairs in a range of indices
  char UpperLower = 'U'; // Use upper triagle of A

  int N   = A.rows();
  double *Amat = new double[N*N];
  int LDA = N;
  double VL = 0.0;
  double VU = 0.0;
  int IL = 1;
  int IU = NumPairs;
  double AbsTolerance = 0.0;
  int NumComputed;
  double *EigVals = new double[N];
  double *EigVecs = new double[N*NumPairs];
  int LDEigVecs = N;
  int *ISuppZ = new int[2*NumPairs];
  int Info;

  // First do workspace query
  int Lwork = -1;
  int LIwork = -1;
  double WorkSize;
  int IWorkSize;
  
  
   F77_DSYEVR (&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
	       &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
	       &LDEigVecs, ISuppZ, &WorkSize, &Lwork, &IWorkSize, &LIwork, 
	       &Info);

   // Now allocate WorkSpace;
   Lwork = (int) floor(WorkSize+0.5);
   LIwork = IWorkSize;
   double *WorkSpace = new double[Lwork];
   int *IWorkSpace = new int[LIwork];
   
   // Copy A int Amat
   for (int row=0; row<N; row++)
     for (int col=0; col<N; col++)
       *(Amat+(col*N)+row) = A(row,col);
  
   F77_DSYEVR (&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
	       &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
	       &LDEigVecs, ISuppZ, WorkSpace, &Lwork, IWorkSpace, &LIwork, 
	       &Info);

   if (Info !=0) 
     {
       fprintf (stderr, "Lapack error in dsyevr_.  Exitting.\n");
       exit(-1);
     }

   // Now copy over output of Vectors and Vals
   Vals.resize(NumPairs);
   Vectors.resize(NumPairs,N);
   
   for (int i=0; i<NumPairs; i++)
     {
       Vals(i) = *(EigVals+i);
       for (int j=0; j<N; j++)
	 Vectors(i,j) = *(EigVecs+(i*N)+j);
     }

   // Now free allocate memory
   delete[] Amat, EigVals, EigVecs, WorkSpace, IWorkSpace, ISuppZ;
   

}




void SymmEigenPairs (const Array<complex<double>,2> &A, int NumPairs,
		     Array<scalar,1> &Vals,
		     Array<complex<double>,2> &Vectors)
{
  char JobType = 'V';    // Find eigenvectors and eignevalues
  char Range   = 'I';    // Find eigenpairs in a range of indices
  char UpperLower = 'U'; // Use upper triagle of A

  int N   = A.rows();
  complex<double> *Amat = new complex<double>[N*N];
  int LDA = N;
  double VL = 0.0;
  double VU = 0.0;
  int IL = 1;
  int IU = NumPairs;
  double AbsTolerance = 0.0;
  int NumComputed;
  double *EigVals = new double[N];
  complex<double> *EigVecs = new complex<double>[N*NumPairs];
  int LDEigVecs = N;
  int *ISuppZ = new int[2*NumPairs];
  int Info;

  // First do workspace query
  int Lwork = -1;
  int LIwork = -1;
  int LRwork = -1;
  complex<double> WorkSize;
  double RWorkSize;
  int IWorkSize;
  
  
   F77_ZHEEVR(&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
	      &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
	      &LDEigVecs, ISuppZ, &WorkSize, &Lwork, &RWorkSize, &LRwork, 
	      &IWorkSize, &LIwork, &Info);
//    fprintf (stderr, "WorkSize  = %1.8f\n", WorkSize.real());
//    fprintf (stderr, "RWorkSize = %1.8f\n", RWorkSize);
//    fprintf (stderr, "IWorkSize = %d\n", IWorkSize);

   // Now allocate WorkSpace;
   Lwork = (int) floor(WorkSize.real()+0.5);
   LIwork = IWorkSize;
   LRwork = (int) floor (RWorkSize+0.5);
   complex<double> *WorkSpace = new complex<double>[Lwork];
   double * RWorkSpace = new double[LRwork];
   int *IWorkSpace = new int[LIwork];
   
   // Copy A int Amat
   for (int row=0; row<N; row++)
     for (int col=0; col<N; col++)
       *(Amat+(col*N)+row) = A(row,col);
   
   F77_ZHEEVR(&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
	      &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
	      &LDEigVecs, ISuppZ, WorkSpace, &Lwork, RWorkSpace, &LRwork, 
	      IWorkSpace, &LIwork, &Info);

   if (Info !=0) {
     fprintf (stderr, "Lapack error in zheevr.  Exitting.\n");
     exit(-1);
   }

   // Now copy over output of Vectors and Vals
   Vals.resize(NumPairs);
   Vectors.resize(NumPairs,N);
   
   for (int i=0; i<NumPairs; i++) {
     Vals(i) = *(EigVals+i);
     for (int j=0; j<N; j++)
       Vectors(i,j) = *(EigVecs+(i*N)+j);
   }

   // Now free allocate memory
   delete[] Amat, EigVals, EigVecs, WorkSpace, IWorkSpace, ISuppZ;
   

}


/// This function returns the determinant of A and replaces A with its
/// cofactors.
double 
DetCofactors (Array<double,2> &A, Array<double,1> &work)
{
  const int maxN = 2000;
  int ipiv[maxN];
  int N = A.rows();
  int M = A.cols();
  assert (N == M);
  assert (N <= maxN);
  // First, transpose A for fortran ordering
//   for (int i=0; i<N; i++)
//     for (int j=0; j<i; j++) {
//       double tmp = A(i,j);
//       A(i,j) = A(j,i);
//       A(j,i) = tmp;
//     }
  Transpose(A);
  
  int info;
  // Do LU factorization
  F77_DGETRF (&N, &M, A.data(), &N, ipiv, &info);
  double det = 1.0;
  int numPerm = 0;
  for (int i=0; i<N; i++) {
    det *= A(i,i);
    numPerm += (ipiv[i] != (i+1));
  }
  if (numPerm & 1)
    det *= -1.0;
  
  int lwork = work.size();
  // Now, do inverse
  F77_DGETRI (&N, A.data(), &N, ipiv, work.data(), &lwork, &info);

  // Now, we have the transpose of Ainv.  Now, just multiply by det:
  A = det * A;
  // And we're done!
  return det;
}

int 
DetCofactorsWorksize(int N)
{
  double work;
  double dummy;
  int info;
  int ipiv;
  int lwork = -1;
  
  F77_DGETRI(&N, &dummy, &N, &ipiv, &work, &lwork, &info);

  return ((int)ceil(work));
}



/// This function returns the determinant of A and replaces A with its
/// cofactors.
complex<double>
ComplexDetCofactors (Array<complex<double>,2> &A, 
		     Array<complex<double>,1> &work)
{
  const int maxN = 2000;
  int ipiv[maxN];
  int N = A.rows();
  int M = A.cols();
  assert (N == M);
  assert (N <= maxN);
  // First, transpose A for fortran ordering
//   for (int i=0; i<N; i++)
//     for (int j=0; j<i; j++) {
//       double tmp = A(i,j);
//       A(i,j) = A(j,i);
//       A(j,i) = tmp;
//     }
  Transpose(A);
  
  int info;
  // Do LU factorization
  F77_ZGETRF (&N, &M, A.data(), &N, ipiv, &info);
  complex<double> det = 1.0;
  int numPerm = 0;
  for (int i=0; i<N; i++) {
    det *= A(i,i);
    numPerm += (ipiv[i] != (i+1));
  }
  if (numPerm & 1)
    det = -det;
  
  int lwork = work.size();
  // Now, do inverse
  F77_ZGETRI (&N, A.data(), &N, ipiv, work.data(), &lwork, &info);

  // Now, we have the transpose of Ainv.  Now, just multiply by det:
  A = det * A;
  // And we're done!
  return det;
}

int 
ComplexDetCofactorsWorksize(int N)
{
  complex<double> work;
  complex<double> dummy;
  int info;
  int ipiv;
  int lwork = -1;
  
  F77_ZGETRI(&N, &dummy, &N, &ipiv, &work, &lwork, &info);

  return ((int)ceil(work.real()));
}
