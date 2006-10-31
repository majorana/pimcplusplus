#include "CubicBspline.h"
#include <iostream>
#include <cstdlib>

void
CubicBspline::Set(double start, double end, Array<double,1> &data,
		  bool interpolating, bool periodic)
{
  Periodic = periodic;
  Interpolating = interpolating;
  
  if (!periodic) {
    cerr << "Non-periodic B-splines not yet implemented.\n";
    abort();
  }
  if (interpolating) {
    cerr << "Interpolating B-splines not yet implemented.\n";
    abort();
  }
  
  P.resize(data.size());
  P = data;
  GridStart = start;
  GridEnd   = end;
  GridDelta = (end-start)/(double)(data.size());
  GridDeltaInv = 1.0/GridDelta;
}


CubicBspline::CubicBspline()
{
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;
}
