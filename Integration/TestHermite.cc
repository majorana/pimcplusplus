#include "HermiteQuad.h"

class TestFunc
{
public:
  inline double operator()(double x)
  { 
    //return ((x+1)*(x+1)*exp(-(2.0*x+1.0)));
    return (x*x);
  }
};


void TestHermite()
{
  HermiteQuadClass<Hermite30, TestFunc> Hermite;
  TestFunc Func;

  double a = 1.0;
  a = 1.0;
  double sigma = sqrt(1.0/(a+a));
  Hermite.SetSigma(sigma);

  double val;
  for (int i=0; i<1000; i++)
    val = Hermite.Integrate(Func);
  double exact = 0.5 * sqrt(M_PI/(a*a*a));
  fprintf (stderr, "Integration = %1.16e\n", val);
  fprintf (stderr, "Exact       = %1.16e\n", exact);
}

class TestFunc3D
{
public:
  inline double operator()(double x, double y, double z)
  {
    //return exp((x*x+y*y+z*z));
    return 1.0;
  }
};

#include <time.h>
void TestHermite3D()
{
  Hermite3DQuadClass<Hermite20, TestFunc3D> Hermite;
  TestFunc3D Func;

  double a = 1.0;
  a = 2.0;
  double sigma = sqrt(1.0/(2.0*a));
  a = 1.0;
  Hermite.SetSigma(sigma);
  double val;
  clock_t Start = clock();
  for (int i=0; i<100000; i++)
    val = Hermite.Integrate(Func);
  clock_t End = clock();
  double time = (double)(End-Start)/1.0e11;
  fprintf (stderr, "Time per integration = %1.4e\n", time);
  double exact = sqrt(M_PI*M_PI*M_PI/(a*a*a));
  fprintf (stderr, "Integration = %1.16e\n", val);
  fprintf (stderr, "Exact       = %1.16e\n", exact);
}

main()
{
  TestHermite();
  TestHermite3D();
}

