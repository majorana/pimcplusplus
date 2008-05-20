#include "BirchEOS.h"

int main()
{
  BirchEOSClass<4> birch;
  TinyVector<double, 6> params(-50.0, 40.0, 3.0, 2000.1, 5, 0.8);
  TinyVector<double, 6> grad, gradFD;


  birch.SetParams(params);
  
  double P1 = birch.Pressure(41.0);
  double P2 = birch.PressureFD(41.0);

  cerr << "P1 = " << P1 << endl;
  cerr << "P2 = " << P2 << endl;

  grad = birch.Grad(41.0);
  gradFD = birch.GradFD(41.0);
  cerr << "grad   = " << grad << endl;
  cerr << "gradFD = " << gradFD << endl << endl;

  double B0 = birch.GetB0();
  double B  = birch.GetB(params[0]);
  double B0FD= birch.GetB0FD();

  cerr << "B0(V0) = " << B << endl;
  cerr << "B0     = " << B0   << endl;
  cerr << "B0FD   = " << B0FD << endl << endl;

  double B0p = birch.GetB0p();
  double B0pFD = birch.GetB0p();
  
  cerr << "B0p   = " << B0p << endl;
  cerr << "B0pFD = " << B0p << endl;

}
