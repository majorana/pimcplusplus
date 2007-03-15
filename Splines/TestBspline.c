#include "bspline.h"
#include <stdio.h>
#include <stdlib.h>

void
Test_1d_s()
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end   = 3.0;
  grid.num = 11;
  float data[] = { 3.0, -4.0, 2.0, 1.0, -2.0, 0.0, 3.0, 2.0, 0.5, 1.0, 3.0 };
  BCtype_s bc;
  bc.lCode = DERIV2; bc.lVal = 10.0;
  bc.rCode = DERIV2; bc.rVal = -10.0;
  
  UBspline_1d_s *spline = (UBspline_1d_s*) create_UBspline_1d_s (grid, bc, data);
  for (double x=1.0; x<=3.00001; x+=0.001) {
    float val, grad, lapl;
    eval_UBspline_1d_s_vgl (spline, x, &val, &grad, &lapl);
    fprintf (stdout, "%1.5f %20.14f %20.14f %20.14f\n", x, val, grad, lapl);
  }

}


void
Test_2d_s()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 10;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 10;
  
  float *data = malloc (x_grid.num * y_grid.num * sizeof(float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = -1.0 + 2.0*drand48();
  BCtype_s x_bc, y_bc;
  x_bc.lCode = PERIODIC; x_bc.lVal = 10.0;
  x_bc.rCode = PERIODIC; x_bc.rVal = -10.0;
  y_bc.lCode = PERIODIC; y_bc.lVal = 10.0;
  y_bc.rCode = PERIODIC; y_bc.rVal = -10.0;
  
  UBspline_2d_s *spline = (UBspline_2d_s*) create_UBspline_2d_s (x_grid, y_grid, x_bc, y_bc, data); 

//   for (int ix=0; ix<x_grid.num+3; ix++) {
//     for (int iy=0; iy<y_grid.num+3; iy++)
//       fprintf (stdout, "%20.14f ", spline->coefs[ix*spline->x_stride + iy]);;
//     fprintf (stdout, "\n");
//   }
  for (double x=x_grid.start; x<=x_grid.end; x+=0.02) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.02) {
      float val;
      eval_UBspline_2d_s (spline, x, y, &val);
      fprintf (stdout, "%20.14f ", val);
    }
    fprintf (stdout, "\n");
  }
}


main()
{
  //  Test_1d_s();
  Test_2d_s();
}
