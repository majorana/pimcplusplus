#include <math.h>
#include <stdio.h>

main()
{
  double box[3] = { 2.0, 2.0, 2.0 };
  double kPrim[3];
  kPrim[0] = 2.0*M_PI/box[0];
  kPrim[1] = 2.0*M_PI/box[1];
  kPrim[2] = 2.0*M_PI/box[2];
  const int iMax = 10;
  
  double beta = 1.0;
  double lambda = 0.5;

  double k[3];
  double Esum = 0.0, Zsum = 0.0;
  for (int ix=-iMax; ix<=iMax; ix++) {
    k[0] = kPrim[0] * ix;
    double Ex = k[0]*k[0]*lambda;
    for (int iy=-iMax; iy<iMax; iy++) { 
      k[1] = kPrim[1] * iy;
      double Ey = k[1]*k[1]*lambda;
      for (int iz=-iMax; iz<iMax; iz++) { 
	k[2] = kPrim[2] * iz;
	double Ez = k[2]*k[2]*lambda;
	double E = Ex + Ey + Ez;
	Esum += E * exp(-beta * E);
	Zsum += exp(-beta*E);
      }
    }
  }
  fprintf (stderr, "Energy = %1.12e\n", Esum/Zsum);
}
