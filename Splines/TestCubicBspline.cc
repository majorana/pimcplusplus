#include "CubicBspline.h"

void TestValue()
{
  Array<double,1> d(10);
  for (int i=0; i<d.size(); i++) {
    double x = 2.0*M_PI*i/(d.size());
    d(i) = sin(x);
  }
  
  CubicBspline spline;
  spline.Set (0.0, 2.0*M_PI, d, false, true);

  for (double x=0.0; x<2.0*M_PI; x+=0.0001) 
    fprintf (stdout, "%20.16e %20.16e\n", x, spline(x));
}

main() 
{
  TestValue();

}
    
    
