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


void LinFit (Array<double,1> &y, Array<double,1> &sigma,    // inputs
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
