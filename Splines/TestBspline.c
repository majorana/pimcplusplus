#include "bspline.h"
#include <stdio.h>

void
Test_1d_s()
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end   = 3.0;
  grid.num = 11;
  float data[] = { 3.0, -4.0, 2.0, 1.0, -2.0, 0.0, 3.0, 2.0, 0.5, 1.0, 3.0 };
  BCtype_s bc;
  bc.lCode = DERIV1; bc.lVal = 10.0;
  bc.rCode = DERIV1; bc.rVal = -10.0;
  
  UBspline_1d_s *spline = (UBspline_1d_s*) create_UBspline_1d_s (grid, bc, data);
  for (double x=1.0; x<=3.00001; x+=0.001) {
    float val;
    eval_UBspline_1d_s (spline, x, &val);
    fprintf (stdout, "%1.5f %20.14f\n", x, val);
  }

}


main()
{
  Test_1d_s();

}
