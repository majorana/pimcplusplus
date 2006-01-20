#include "FastCubicSpline.h"
#include "CubicSpline.h"
#include "MyTricubicSpline.h"

#include<time.h>

void
FastCubicSpeed()
{
  Array<double,1> y(500);
  for (int i=0; i<y.size(); i++) {
    double x = 2.0*M_PI*(double)i/(double)(y.size()-1);
    y(i) = sin(x);
  }
  FastCubicSpline spline;
  spline.Init (0.0, 2.0*M_PI, y, 1.0, 1.0);

  clock_t start, end, rstart, rend;

  const int numSamples = 100000000;  
  rstart = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
  }
  rend = clock();

  start = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double p = spline(x);
  }
  end = clock();

  double time = (double)((end-start)-(rend-rstart))/(double)CLOCKS_PER_SEC;


  fprintf (stderr, "Time for 1 spline evaluation = %1.5e\n",
	   time/(double)numSamples);
}


void
CubicSpeed()
{
  Array<double,1> y(500);
  LinearGrid xGrid (0.0, 2.0*M_PI, 500);
  for (int i=0; i<y.size(); i++) {
    double x = 2.0*M_PI*(double)i/(double)(y.size()-1);
    y(i) = sin(x);
  }
  CubicSpline spline (&xGrid, y);

  clock_t start, end, rstart, rend;

  const int numSamples = 100000000;  
  rstart = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
  }
  rend = clock();

  start = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double p = spline(x);
  }
  end = clock();

  double time = (double)((end-start)-(rend-rstart))/(double)CLOCKS_PER_SEC;


  fprintf (stderr, "Time for 1 spline evaluation = %1.5e\n",
	   time/(double)numSamples);
}

void
TricubicSpeed()
{
  Array<double,1> y(500);
  LinearGrid xGrid (0.0, 2.0*M_PI, 500);
  LinearGrid yGrid (0.0, 1.0, 50);
  LinearGrid zGrid (0.0, 1.0, 50);
  Array<double,3> F(500,50,50);
  for (int i=0; i<F.extent(0); i++) {
    double x = xGrid(i);
    for (int j=0; j<F.extent(1); j++) {
      double y = yGrid(j);
      for (int k=0; k<F.extent(2); k++) {
	double z = zGrid(k);
	F(i,j,k) = sin(x)*(1.0+y*y*z*z*z);
      }
    }
  }
  MyTricubicSpline spline (&xGrid, &yGrid, &zGrid, F);

  clock_t start, end, rstart, rend;

  const int numSamples = 10000000;  
  rstart = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double y = drand48();
    double z = drand48();
  }
  rend = clock();

  start = clock();
  for (int i=0; i<numSamples; i++) {
    double x = 2.0*M_PI*drand48();
    double y = drand48();
    double z = drand48();
    double p = spline(x,y,z);
  }
  end = clock();

  double time = (double)((end-start)-(rend-rstart))/(double)CLOCKS_PER_SEC;


  fprintf (stderr, "Time for 1 spline evaluation = %1.5e\n",
	   time/(double)numSamples);
}


main()
{
  FastCubicSpeed();
  CubicSpeed();
  TricubicSpeed();
}
