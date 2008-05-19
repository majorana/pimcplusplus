#include "BirchEOS.h"

int main()
{
  BirchEOSClass<2> birch;
  TinyVector<double, 4> params(-50.0, 40.0, 3.0, 2.1);
  TinyVector<double, 4> grad, gradFD;


  birch.SetParams(params);
  
  double P1 = birch.Pressure(41.0);
  double P2 = birch.PressureFD(41.0);

  cerr << "P1 = " << P1 << endl;
  cerr << "P2 = " << P2 << endl;

  grad = birch.Grad(41.0);
  gradFD = birch.GradFD(41.0);
  cerr << "grad   = " << grad << endl;
  cerr << "gradFD = " << gradFD << endl;


}
