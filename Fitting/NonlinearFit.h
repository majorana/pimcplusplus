#ifndef NONLINEAR_FIT_H
#define NONLINEAR_FIT_H

#include <blitz/tinyvec.h>
#include <blitz/array.h>

using namespace blitz;

// Template class for fitting a function with M nonlinear parameters.
template<int M, typename ModelType>
class NonlinearFitClass
{
private:
  ModelType &Model;
  TinyVector<double,M> Params, Beta;
  TinyMatrix<double,M,M> Alpha;
  double Chi2 (const Array<double,1> &x, const Array<double,1> &y,
	       const Array<double,1> &y, TinyVector<double,M> params);
  void CalcAlphaBeta(const Array<double,1> &x, 
		     const Array<double,1> &y,
		     const Array<double,1> &sigma,
		     TinyVector<double,M> params);
public:
  void Fit(const Array<double,1> &x, const Array<double,1> &y,
	   const Array<double,1> &sigma, TinyVector<double,M> &params);

  NonlinearFitClass(ModelType model) :
    Model(model)
  {
    
  };
};

#endif
