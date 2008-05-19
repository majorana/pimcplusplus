#include "BirchEOS.h"

int main()
{
  BirchEOSClass<2> birch;
  TinyVector<double, 4> params(-50.0, 40.0, 3.0, 2.1);

  birch.SetParams(params);
  
  double P1 = birch.Pressure(41.0);
  double P2 = birch.PressureFD(41.0);

  cerr << "P1 = " << P1 << endl;
  cerr << "P2 = " << P2 << endl;


}
