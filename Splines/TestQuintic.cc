#include "QuinticSpline.h"

void TestQuintic1()
{
  Array<double,1> Y(6);
  LinearGrid X(0.0, 5.0*M_PI, 6);
  
  Y(0)= 0.0; Y(1)=1; Y(2)=0.0; Y(3)=1.0; Y(4) = 0.0; Y(5) = 1.0;
  double startDeriv = 0.0;
  double endDeriv = 0.0;
  double startDeriv2 = 0.0;
  double endDeriv2 = -0.0;
  
  QuinticSpline sp(&X, Y, startDeriv, endDeriv, startDeriv2, endDeriv2);

  FILE *fout = fopen ("quintic1.dat", "w");
  LinearGrid X2(0.0, 5.0*M_PI, 1000);
  for (int i=0; i<X2.NumPoints; i++) {
    double x = X2(i);
    double y = sp(x);
    double dy = sp.Deriv(x);
    double d2y = sp.Deriv2(x);
    double d3y = sp.Deriv3(x);
    double d4y = sp.Deriv4(x);
    fprintf (fout, "%1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", 
	     x, y, dy, d2y, d3y, d4y);
  }
  fclose (fout);

}


void TestQuintic2()
{
  Array<double,1> Y(6);
  LinearGrid X(0.0, 5.0, 6);
  
  Y(0)= 0.0; Y(1)=1.0; Y(2)=0.0; Y(3)=1.0; Y(4) = 0.0; Y(5) = 1.0;
  
  QuinticSpline sp(&X, Y);

  FILE *fout = fopen ("quintic2.dat", "w");
  LinearGrid X2(0.0, 5.0, 1000);
  for (int i=0; i<X2.NumPoints; i++) {
    double x = X2(i);
    double y = sp(x);
    double dy = sp.Deriv(x);
    double d2y = sp.Deriv2(x);
    double d3y = sp.Deriv3(x);
    double d4y = sp.Deriv4(x);
    fprintf (fout, "%1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", 
	     x, y, dy, d2y, d3y, d4y);
  }
  fclose (fout);

}



main()
{
  TestQuintic1();
  TestQuintic2();
}
