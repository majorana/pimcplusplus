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

class NormIntegrand
{
  WeightFuncClass &WF;
  PolynomialClass &P;
  double xMin;
public:

  inline double operator()(double delta)
  {
    double x = xMin + delta;
    double p = P(x);
    return (p*p*WF(x));
  }
  NormIntegrand(double xmin, PolynomialClass &p,
		WeightFuncClass &wf) : xMin(xmin), WF(wf), P(p) 
  {} 
};



void PolynomialSetClass::MakeOrthoSet(int order, double xmin, double xmax,
				      WeightFuncClass &wf)
{
  const double absTolerance = 1.0e-8*(xmax-xmin);
  const double relTolerance = 1.0e-7;
  /// First, calculate S_n = \int_{xmin}^{xmax} x^n w(x) \ dx
  Array<double,1> Sn(2*order+1);
  for (int n=0; n<=(2*order); n++) {
    PolyIntegrand polyInt(n, wf);
    GKIntegration<PolyIntegrand,GK15> Integrator(polyInt);
    // Integrate until the absolute tolerance OR the relative
    // tolerance is satisfied.
    Sn(n) = Integrator.Integrate(xmin, xmax, 
				 absTolerance, relTolerance, false);
  }
    
  P.resize(order+1);
  for (int i=0; i<=order; i++) {
    P(i).SetOrder(order);
    for (int j=0; j<=order; j++)
      P(i)[j] = 0.0;
  }

  // Initialize first function to a constant
  P(0)[0] = 1.0;
  PolynomialClass P2 = P(0)*P(0);
  NormIntegrand normInt(xmin, P(0), wf);
  GKIntegration<NormIntegrand,GK15> Integrator(normInt);
  Integrator.SetRelativeErrorMode();
  double norm = Integrator.Integrate(0.0, xmax-xmin, relTolerance);
  P(0) = (1.0/sqrt(norm))*P(0);

//   double norm=0.0;
//   for (int i=0; i<=(2*order); i++)
//     norm += P2[i]*Sn(i);
//   P(0) = (1.0/sqrt(norm))*P(0);


  // Now construct remaining functions
  for (int n=1; n<=order; n++) {
    // First, set P(n) = x*P(n-1)
    P(n)[0] = 0.0;
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

    for (int j=0; j<=order; j++)
      if (isnan(P(n)[j])) {
	cerr << "Error in PolynomialSet construction.\n";
	cerr << "norm = " << norm << endl;
      }

    // New normalization:
    NormIntegrand normInt(xmin, P(n), wf);
    GKIntegration<NormIntegrand,GK15> Integrator(normInt);
    //cerr << "<Normalizing>\n";
    Integrator.SetRelativeErrorMode();
    double norm = Integrator.Integrate(0.0, xmax-xmin, relTolerance);
    //cerr << "</Normalizing>\n";
    P(n) = (1.0/sqrt(norm))*P(n);
    for (int j=0; j<=order; j++)
      if (isnan(P(n)[j])) {
	cerr << "Error in PolynomialSet construction.\n";
	cerr << "norm = " << norm << endl;
      }
    // Faster, older normalization -- might give NAN's
//     PolynomialClass P2 = P(n)*P(n);
//     double norm=0.0;
//     for (int i=0; i<=(2*order); i++)
//       norm += P2[i]*Sn(i);
//     P(n) = (1.0/sqrt(norm))*P(n);
//     for (int j=0; j<=order; j++)
//       if (isnan(P(n)[j])) {
// 	cerr << "Error in PolynomialSet construction.\n";
// 	cerr << "norm = " << norm << endl;
//       }

  }
}



