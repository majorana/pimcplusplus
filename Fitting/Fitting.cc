#include "Fitting.h"
#include "../MatrixOps/MatrixOps.h"

inline Array<double,1> operator*(Array<double,2> &A, Array<double,1> &x)
{
  assert (A.cols() == x.size());
  Array<double,1> Ax(A.rows());
  Ax = 0.0;
  for (int i=0; i<A.rows();i++)
    for (int j=0; j<A.cols(); j++)
      Ax(i) += A(i,j) * x(j);
  return Ax;
}

inline Array<double,1> diag (Array<double,2> &A)
{
  assert (A.rows() == A.cols());
  Array<double,1> d(A.rows());
  for (int i=0; i<A.rows(); i++)
    d(i) = A(i,i);
  return (d);
}


void LinFitLU (Array<double,1> &y, Array<double,1> &sigma,    // inputs
	       Array<double,2> &F,                     // input
	       Array<double,1> &a, Array<double,1> &errors) // outputs
{
  int N = F.rows();
  int M = F.cols();

  if (y.rows() != F.rows()) {
    cerr << "Differernt number of rows in basis functions that in data "
	 << "points in LinFit.  Exitting.\n";
    abort();
  }

  assert (y.rows() == sigma.rows());
  a.resize (F.cols());
  errors.resize(F.cols());

  // First, construct A matrix
  Array<double,2> A(F.rows(), F.cols());
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)
      A(i,j) = F(i,j) / sigma(i);
  
  // Next, construct b vector
  Array<double,1> b(y.rows());
  for (int i=0; i<N; i++)
    b(i) = y(i)/sigma(i);
  
  // Next, construct alpha matrix
  Array<double,2> alpha(M,M);
  alpha = 0.0;
  for (int j=0; j<M; j++)
    for (int k=0; k<M; k++)
      for (int i=0; i<N; i++)
	alpha(k,j) += A(i,j) * A(i,k);
  
  // Next, construct beta vector
  Array<double,1> beta(M);
  beta = 0.0;
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta(k) += b(i)*A(i,k);

  // Now, invert alpha
  Array<double,2> C(M,M);
  C = Inverse(alpha);
  a = C * beta;
  for (int i=0; i<M; i++)
    errors(i) = sqrt(C(i,i));
  //  errors = sqrt(diag(C));  
}




void LinFitSVD (Array<double,1> &y, Array<double,1> &sigma,  // inputs
		Array<double,2> &F,                          // input
		Array<double,1> &a, Array<double,1> &errors, // outputs
		double tolerance)
{
  int N = F.rows();
  int M = F.cols();

  if (y.rows() != F.rows()) {
    cerr << "Differernt number of rows in basis functions that in data "
	 << "points in LinFit.  Exitting.\n";
    abort();
  }

  assert (y.rows() == sigma.rows());
  a.resize (F.cols());
  errors.resize(F.cols());

  // First, construct A matrix
  Array<double,2> A(F.rows(), F.cols());
  for (int i=0; i<N; i++)
    for (int j=0; j<M; j++)
      A(i,j) = F(i,j) / sigma(i);
  
  // Next, construct b vector
  Array<double,1> b(y.rows());
  for (int i=0; i<N; i++)
    b(i) = y(i)/sigma(i);
  
  // Now, compute SVD
  Array<double,2> U,V;
  Array<double,1> S;
  SVdecomp (A, U, S, V);

  //  cerr << "S = " << S << endl;
  // Zero out near-singular values
  double Smax=S(0);
  for (int i=1; i<S.size(); i++)
    Smax = max (S(i),Smax);
  Array<double,1> Sinv(S.size());
  for (int i=0; i<S.size(); i++)
    Sinv(i) = (S(i) < (tolerance*Smax)) ? 0.0 : (1.0/S(i));
  
  a = 0.0;
  for (int k=0; k<a.size(); k++) 
    for (int i=0; i<U.cols(); i++) {
      double coef = 0.0;
      for (int j=0; j < U.rows(); j++)
	coef += U(j,i) * b(j);
      coef *= Sinv(i);
      a(k) += coef * V(k,i);
    }
  
  errors = 0.0;
  for (int j=0; j<errors.size(); j++) {
    for (int i=0; i<V.cols(); i++)
      errors(j) += V(j,i)*V(j,i)*(Sinv(i)*Sinv(i));
    errors(j) = sqrt(errors(j));
  }
}
