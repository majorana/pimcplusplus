#include "LegendrePoly.h"
#include "SpecialFunctions.h"


void LegendreTest()
{
  Array<double,1> Pn(50);
  for (int n=0; n<50; n++) {
    for (double x=-1.0; x<1.0; x+=0.001) {
      double l1 = Legendre(n,x);
      LegendrePoly(x, Pn);
      double l2 = Pn(n);
      fprintf (stderr, "%7.4f %21.14e %21.14e\n", x, l1, l2);
    }
  }
}

main()
{
  LegendreTest();

}


