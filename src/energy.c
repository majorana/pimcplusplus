#include <math.h>
#include <stdio.h>

#include <iostream>
main()
{
  double L=3.0;
  double lambda = 1.0;
  double Esum, Z;
  Esum = 0.0;
  Z = 0.0;
  double beta=1.0;
  int cut = 3;

  for (int i=-cut; i<=cut; i++) {
    double kx = 2.0*M_PI/L * (double)i;
    for (int j=-cut; j<=cut; j++) {
      double ky = 2.0*M_PI/L * (double)j;
      for (int k=-cut; k<=cut; k++) {
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
  
  double Esum2, Z2;
  Esum2 = 0.0; Z2 = 0.0;
  int numPerDim = 2*cut+1;
  int num = numPerDim*numPerDim*numPerDim;
  for (int i=0; i < num; i++) {
    int ix = i/(numPerDim*numPerDim) + cut;
    int iy = (i/numPerDim)%numPerDim + cut;
    int iz = i % numPerDim + cut;
    double kx = 2.0*M_PI/L * ix;
    double ky = 2.0*M_PI/L * iy;
    double kz = 2.0*M_PI/L * iz;
    cerr << ix << ", " << iy << ", " << iz <<endl;
    double E = lambda*(kx*kx+ky*ky+kz*kz);
    Esum2 += exp(-beta*E)*E;
    Z2 += exp(-beta*E)*E;
  }
	
  fprintf (stderr, "Eavg = %1.8f\n", Esum/Z);
  fprintf (stderr, "Eavg = %1.8f\n", Esum2/Z2);
}
