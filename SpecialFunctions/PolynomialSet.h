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
  inline double operator()(int n, double x)
  { return P(n)(x);}
  inline PolynomialClass& operator()(int n)
  { return P(n);}

  void MakeOrthoSet(int order, double xmin, double xmax, 
		    WeightFuncClass &wf);
};


#endif
