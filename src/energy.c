#include <math.h>
#include <stdio.h>

main()
{
  double L=3.0;
  double lambda = 1.0;
  double Esum, Z;
  Esum = 0.0;
  Z = 0.0;
  double beta=2.0;

  for (int i=-200; i<=200; i++) {
    double kx = 2.0*M_PI/L * (double)i;
    for (int j=-200; j<=200; j++) {
      double ky = 2.0*M_PI/L * (double)j;
      for (int k=-200; k<=200; k++) {
	double kz = 2.0*M_PI/L * (double) k;
	double E = lambda*(kx*kx+ky*ky+kz*kz);
	//	fprintf (stderr, "E = %1.10f\n", E);
	//if (kx!=0 || ky!=0 || kz!=0) {
	  Esum += exp(-beta*E)*E;
	  Z += exp(-beta*E);
	  //	}
      }
    }
  }
	
  fprintf (stderr, "Eavg = %1.8f\n", Esum/Z);
}
