#ifndef HERMITE_QUAD_H
#define HERMITE_QUAD_H

#include "math.h"

class Hermite20
{
public:
  static const int n=20;
  static const double u[20];
  static const double weight[20];
  double w[20];
  double x[20];
};

class Hermite30
{
public:
  static const int n=30;
  static const double u[30];
  static const double weight[30];
  double w[30];
  double x[30];
};



template <class RuleClass, class ScaledIntegrand> 
class HermiteQuadClass
{
 private:
  RuleClass Rule;
  double sigma;
 public:
  void SetSigma(double mysigma)
    {
      sigma = mysigma;
      for (int i=0; i<Rule.n; i++)
	Rule.x[i] = Rule.u[i]*M_SQRT2*sigma;
    }  

  inline double Integrate(ScaledIntegrand &Integrand)
    {
      double sum=0; 
      for (int i=0; i<Rule.n; i++)
	sum += Rule.w[i]*Integrand(Rule.x[i]);
      sum *= M_SQRT2*sigma;
      return (sum);
    }
  HermiteQuadClass()
  {
    for (int i=0; i<Rule.n; i++)
      Rule.w[i] = Rule.weight[i]*exp(-Rule.u[i]*Rule.u[i]);
  }
};


void TestHermite();

#endif
