#include "Polynomial.h"

void TestPoly()
{
  PolynomialClass a(3), b(2), c(4);
  a[0]= 2.0; a[1]= -1.0; a[2] = 0.2; a[3]=1.5;
  b[0]= 0.5; b[1]= 2.1;  b[2]=-4.9;
  c[0]= 0.3; c[1]= 1.5;  c[2]=1.5; c[3]=-2.0; c[4]=0.2;
 
  double x = 3.14159;
  double aofx = a(x);
  double bofx = b(x);
  double cofx = c(x);
  
  PolynomialClass ab(a*b);
  PolynomialClass aplusb;
  aplusb = a+b;
 
  fprintf (stderr, "a(x)*b(x) = %1.16e\n", a(x)*b(x));
  fprintf (stderr, "(a*b)(x)  = %1.16e\n", ab(x));
  fprintf (stderr, "a(x)+b(x) = %1.16e\n", a(x)+b(x));
  fprintf (stderr, "(a+b)(x)  = %1.16e\n", aplusb(x));

}

main()
{
  TestPoly();
}
