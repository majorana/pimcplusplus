#include "PolynomialSet.h"
#include "../Integration/GKIntegration.h"

class PolyIntegrand
{
  int n;
  WeightFuncClass &WF;
public:

  inline double operator()(double x)
  {
    double x2n = 1.0;
    for (int i=0; i<n; i++)
      x2n *= x;
    return (WF(x) * x2n);
  }
  PolyIntegrand(int n_, WeightFuncClass &wf) : n(n_), WF(wf) 
  {} 
};



void PolynomialSetClass::MakeOrthoSet(int order, double xmin, double xmax,
				      WeightFuncClass &wf)
{
  const double tolerance = 1.0e-9;
  /// First, calculate S_n = \int_{xmin}^{xmax} x^n w(x) \ dx
  Array<double,1> Sn(2*order+1);
  for (int n=0; n<=(2*order); n++) {
    PolyIntegrand polyInt(n, wf);
    GKIntegration<PolyIntegrand,GK15> Integrator(polyInt);
    Integrator.SetRelativeErrorMode();
    Sn(n) = Integrator.Integrate(xmin, xmax, tolerance);
  }
    

  P.resize(order+1);
  for (int i=0; i<=order; i++) {
    P(i).SetOrder(order);
    for (int j=0; j<=order; j++)
      P(i)[j] = 0.0;
  }

  // Initialize first function to a constant
  P(0)[0] = 1.0/sqrt(xmax-xmin);
  for (int n=1; n<=order; n++) {
    // First, set P(n) = x*P(n-1)
    for (int j=1; j<=order; j++)
      P(n)[j] = P(n-1)[j-1];
    // Now, orthogonalize
    for (int m=0; m<n; m++) {
      PolynomialClass Pnm = P(n)*P(m);
      double overlap = 0.0;
      for (int i=0; i<=(2*order); i++)
	overlap += Pnm[i]*Sn(i);
      P(n) = P(n) - overlap*P(m);
    }
    // Now normalize
    PolynomialClass P2 = P(n)*P(n);
    double norm=0.0;
    for (int i=0; i<=(2*order); i++)
      norm += P2[i]*Sn(i);
    P(n) = (1.0/sqrt(norm))*P(n);
  }
}



