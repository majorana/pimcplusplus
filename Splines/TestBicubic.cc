#include "BicubicSpline.h"

void BCSplineTest1()
{
  LinearGrid xgrid(0.0, 4.0*M_PI, 200);
  LinearGrid ygrid(0.0, 4.0*M_PI, 100);
  
  Array<double,2> F(200,100);
  for (int i=0; i<200; i++) {
    double x = xgrid(i);
    for (int j=0; j<100; j++) {
      double y = ygrid(j);
      F(i,j) = sin(x)*sin(y);
    }
  }
  BicubicSpline BCspline;   
  BCspline.Init(&xgrid, &ygrid, F);
  /*FILE *exout = fopen ("BCexact", "w");
  FILE *spout = fopen ("BCspline", "w");
  LinearGrid finegrid(0.0,4.0*M_PI, 1000);
  for (int i=0; i<1000; i++) {
    double x = finegrid(i);
    for (int j=0; j<1000; j++) {
      double y = finegrid(j);      
      double ex = sin(x)*sin(y);
      fprintf (exout, "%1.16e ", ex);
      BCspline(x,y);
      fprintf (spout, "%1.16e ", BCspline(x,y));
    }
    fprintf (exout, "\n");
    fprintf (spout, "\n");
  }
  fclose(exout);
  fclose(spout);*/
  for (int k=0; k<10000000; k++) {
    double x = drand48()*4.0*M_PI;
    double y = drand48()*4.0*M_PI;
    BCspline(x,y);
  }
  


  /*LinearGrid x2(0.0, 10.0*M_PI, 10000);
  for (int i=0; i<10000; i++) {
    double y = x2(i);
    double x = 0.0;
    fprintf (stderr, "%1.17e %1.17e %1.17e %1.17e\n", x,
	     BCspline(0, y), Cspline(y), cos(x) * cos(y));
	     }*/
}

main()
{
  BCSplineTest1();


}
