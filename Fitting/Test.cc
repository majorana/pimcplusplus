#include "Fitting.h"
#include <stdlib.h>
#include <math.h>

void PolyTest()
{
  int N = 100000;
  int M = 4;
  double xmax = 3.0;
  double noise = 0.1;
  Array<double,1> coefs(M);
  Array<double,1> x(N), y(N);
  Array<double,2> F(N,M);

  coefs = 1.0, -2.0, 0.5, 0.2;
  for (int i=0; i<N; i++) {
    x(i) = xmax * (double)i/(double)(N-1);
    y(i) = 0;
    for (int j=0; j<coefs.size(); j++)
      {
	y(i) += coefs(j) * pow(x(i),(double)j);
	F(i,j) = pow(x(i), (double)j);
      }
    y(i) += noise * (-1.0 + 2.0*drand48());
  }
  cerr << "Exact coefs = " << coefs << endl;
  Array<double,1> fitcoefs(M), errors(M), sigma(N);
  sigma = noise;
  LinFit (y, sigma, F, fitcoefs, errors);
  cerr << "Fit Coefs   = " << fitcoefs << endl;
  cerr << "Errors      = " << errors << endl;


}
    
      


main()
{
  PolyTest();
  cerr << "No tests implemented for Fitting.cc.\n";
}
