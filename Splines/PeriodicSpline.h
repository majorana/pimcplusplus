#ifndef PERIODIC_SPLINE_H
#define PERIODIC_SPLINE_H

#include "Grid.h"

class PeriodicSpline
{
private:
  // The first component has y, the second has dy/dx
  Array<Vec2,1> F;
  Grid *grid;
  void Update();
  bool IsUp2Date;
  
public:
  inline double operator()(double x);
  
  void Init (Grid *newGrid, const Array<double,1> &data);
  PeriodicSpline() : IsUp2Date(false)
  {
    // do nothing
  }
};

inline double
PeriodicSpline::operator()(double x)
{
  if (!IsUp2Date)
    Update();
  int i = grid->ReverseMap(x);
  double h = (*grid)(i+1) - (*grid)(i);
  double t = (x - (*grid)(i))/h;
  double tm1 = t-1.0;
  double p1 = (1.0 + 2.0*t)*tm1*tm1;
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;

  return (F(i)[0]*p1 + F(i+1)[0]*p2 + h*(F(i)[1]*q1 + F(i+1)[1]*q2));
}


#endif
