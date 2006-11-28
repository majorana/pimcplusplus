#include "TricubicNUBspline.h"

void
TestLinear()
{
  LinearGrid xGrid (1.0, 5.0, 22);
  LinearGrid yGrid (2.0, 9.0, 31);
  LinearGrid zGrid (3.1, 5.7, 19);

  TricubicNUBspline<double> spline;
  Array<double,3> data (xGrid.NumPoints-1, yGrid.NumPoints-1, zGrid.NumPoints-1);
  for (int ix=0; ix<xGrid.NumPoints-1; ix++)
    for (int iy=0; iy<yGrid.NumPoints-1; iy++)
      for (int iz=0; iz<zGrid.NumPoints-1; iz++) 
	data(ix,iy,iz) = 2.0*(drand48()-0.5);
  spline.Init (&xGrid, &yGrid, &zGrid, data);

  for (double x=1.0; x<5.0; x+=0.001) {
    double y = 2.342;
    double z = 5.2341;
    TinyVector<double,3> r(x,y,z);
    double val = spline (r);
    fprintf (stderr, "%20.16e %20.16e\n", x, val);
  }

}

void
TestNonlinear()
{
  OptimalGrid2 xGrid, yGrid, zGrid;  
  xGrid.Init (1.1, 5.5, 3.0, 22);
  yGrid.Init (3.1, 8.3, 4.0, 25);
  zGrid.Init (2.1, 7.9, 5.0, 26);
  //   LinearGrid xGrid, yGrid, zGrid;
  //   xGrid.Init (1.1, 5.5, 22);
  //   yGrid.Init (3.1, 8.3, 25);
  //   zGrid.Init (2.1, 7.9, 26);
//   Array<double,1> points(11);
//   points = 0.0, 1.0, 1.4, 1.8, 6.0, 7.0, 9.0, 10.0, 12.8, 13.0, 15.0;
//   GeneralGrid xGrid, yGrid, zGrid;
//   xGrid.Init (points);
//   yGrid.Init (points);
//   zGrid.Init (points);

//  TricubicNUBspline<double, LinearGrid, LinearGrid, LinearGrid> spline;
//  TricubicNUBspline<double, GeneralGrid, GeneralGrid, GeneralGrid>
//  spline;
  TricubicNUBspline<double, OptimalGrid2, OptimalGrid2, OptimalGrid2> spline;

  Array<double,3> data (xGrid.NumPoints-1, yGrid.NumPoints-1, zGrid.NumPoints-1);
  for (int ix=0; ix<xGrid.NumPoints-1; ix++)
    for (int iy=0; iy<yGrid.NumPoints-1; iy++)
      for (int iz=0; iz<zGrid.NumPoints-1; iz++) 
	data(ix,iy,iz) = 2.0*(drand48()-0.5);
  spline.Init (&xGrid, &yGrid, &zGrid, data, PERIODIC, PERIODIC, PERIODIC);

  for (double x=1.1; x<=5.5; x+=0.001) {
    double y = 3.342;
    double z = 5.2341;
    TinyVector<double,3> r(x,y,z);
    double val = spline (r);
    fprintf (stderr, "%20.16e %20.16e\n", x, val);
  }

}

void
TestComplexNonlinear()
{
  OptimalGrid2 xGrid, yGrid, zGrid;  
  xGrid.Init (1.1, 5.5, 3.0, 22);
  yGrid.Init (3.1, 8.3, 4.0, 25);
  zGrid.Init (2.1, 7.9, 5.0, 26);
  TricubicNUBspline<complex<double>, OptimalGrid2, OptimalGrid2, OptimalGrid2> spline;
  
  Array<complex<double>,3> data (xGrid.NumPoints-1, yGrid.NumPoints-1, zGrid.NumPoints-1);
  for (int ix=0; ix<xGrid.NumPoints-1; ix++)
    for (int iy=0; iy<yGrid.NumPoints-1; iy++)
      for (int iz=0; iz<zGrid.NumPoints-1; iz++) 
	data(ix,iy,iz) = complex<double>(2.0*(drand48()-0.5),2.0*(drand48()-0.5)); 
  spline.Init (&xGrid, &yGrid, &zGrid, data, PERIODIC, PERIODIC, PERIODIC);

  for (double x=1.1; x<=5.5; x+=0.001) {
    double y = 3.342;
    double z = 5.2341;
    TinyVector<double,3> r(x,y,z);
    complex<double> val = spline (r);
    fprintf (stderr, "%20.16e %20.16e %20.16e\n", x, val.real(), val.imag());
  }

}

main()
{
  TestComplexNonlinear();
}
