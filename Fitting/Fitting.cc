#include "Fitting.h"

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
	alpha += A(i,j) * A(i,k);
  
  // Next, construct beta vector
  Array<double,1> beta(M);
  beta = 0.0;
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta += b(i)*A(i,k);

  // Now, invert alpha
  Array<double,2> C(M,M);
  
  C = Inverse(alpha);
}
