#include "PolynomialSet.h"

class ConstClass : public WeightFuncClass
{
public:
  double operator()(double x)
  { return 1.0; }
};

void TestPolySet()
{
  const int order=10;
  PolynomialSetClass Pset;
  ConstClass wf;
  Pset.MakeOrthoSet(order, -1.0, 1.0, wf);
  for (double x=-1.0; x<=1.00001; x+=0.001) {
    fprintf (stderr, "%1.4e ", x);
    for (int n=0; n<=order; n++)
      fprintf (stderr, "%1.14e ", Pset(n, x));
    fprintf (stderr, "\n");
  }
}


main()
{
  TestPolySet();
}
