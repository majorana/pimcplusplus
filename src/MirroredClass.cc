#include "MirroredClass.h"

ModeType ActiveCopy=OLDMODE;

void MirroredClassTest()
{
  MirroredClass<double> a;
  double b;

  b = 1.0;
  a = b;
  cerr << "a = " << a << endl;;

  Mirrored1DClass<double> c;
  Array<double,1> d;
  c.resize(3);
  d.resize(3);
  c(0) = 1.0; 
  c(1) = 2.0; c(2) = 3.0;
  d = c;
  cerr << "c = " << c.data() << endl;


  Mirrored2DClass<double> e;
  Array<double,2> f;

  e.resize(2,3);
  e.data() = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  f = e;

}
