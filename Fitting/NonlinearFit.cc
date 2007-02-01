#include "NonlinearFit.h"

template<int M, typename ModelType> double
NonlinearFitClass<M,ModelType>::Chi2(const Array<double,1> &x, 
				     const Array<double,1> &y,
				     const Array<double,1> &sigma,
				     TinyVector<double,M> params)
{
  int N = x.size();
  assert (y.size()     == N);
  assert (sigma.size() == N);

  double chi2 = 0;
  for (int i=0; i<N; i++) {
    double val = Model (x(i), Params);
    chi2 += (val-y(i))*(val-y(i))/(sigma(i)*sigma(i));
  }
  return chi2;
}

template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::CalcAlphaBeta (const Array<double,1> &x, 
					       const Array<double,1> &y,
					       const Array<double,1> &sigma,
					       TinyVector<double,M> params)
{
  int N = x.size();
  assert (y.size()     == N);
  assert (sigma.size() == N);
  
  for (int k=0; k<M; k++) {
    Beta[k] = 0.0;
    fori (int l=0; l<M; l++)
      Alpha(k,l) = 0.0;
  }
  
  for (int i=0; i<N; i++) {
    double val = Model(x(i), params);
    TinyVector<double,M> grad = Model.Grad(x(i),params);
    for (int k=0; k<M; k++)
      Beta[k] += grad[k]*(y(i)-val)/(sigma(i)*sigma(i));
  }
  for (int i=0; i<N; i++) {
    double val = Model(x(i), params);
    TinyVector<double,M> grad = Model.Grad(x(i),params);
    for (int k=0; k<M; k++)
      for (int l=0; l<M; l++)
	Alpha(k,l) += grad[k]*grad[l]/(sigma(i)*sigma(i));
  }
}


template<int M, typename ModelType> void
NonlinearFitClass<M,ModelType>::Fit (const Array<double,1> &x,
				     const Array<double,1> &y,
				     const Array<double,1> &sigma,
				     TinyVector<double,M> &params)
{
  double chiNow = Chi2 (x, y, sigma, params);
  double lambda = 0.001;
  
  bool done = false;
  while (!done) {
    CalcAlphaBeta (x, y, sigma, params);
    for (int i=0; i<M; i++)
      Alpha(i,i) *= (1.0+lambda);
    
  }


}
