#ifndef POLYNOMIAL_SET_CLASS_H
#define POLYNOMIAL_SET_CLASS_H
#include "Polynomial.h"

class WeightFuncClass 
{
public:
  virtual double operator()(double x) = 0;
};


class PolynomialSetClass
{
  Array<PolynomialClass,1> P;
public:
  inline double operator()(int k, double x)
  { return P(k)(x);}
  inline PolynomialClass& operator()(int k)
  { return P(k);}

  void MakeOrthoSet(int order, double xmin, double xmax, 
		    WeightFuncClass &wf);
};


#endif
