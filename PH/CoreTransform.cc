#include "../Integration/Integrate.h"
#include "CoreTransform.h"

double dx_dr (double r, double x, void *potptr)
{
  Potential &pot = *(Potential *)potptr;
  double A = pot.A(r);

  return (1.0/sqrt(A));
}

double dr_dx (double x, double r, void *potptr)
{
  Potential &pot = *(Potential *)potptr;
  double A = pot.A(r);

  return (sqrt(A));
}

void
CoreTransform::Initialize(Potential *pot, int NumPoints)
{

  // First, we calculate the r->x transform by integrating
  rMax = pot->CoreRadius();
  rgrid.Init (0.0, rMax, NumPoints);
  Array<double,1> Temp(NumPoints);

  Temp(0) = 0.0;
  IntegrateFirstOrderNS(rgrid, 0, NumPoints-1, Temp, dx_dr, pot);
  xMax = Temp(NumPoints-1);

  double A0 = pot->A(0.0);
  double StartDeriv = 1.0/sqrt(A0);
  double EndDeriv = 1.0;

  x_of_r.Init(&rgrid, Temp, StartDeriv, EndDeriv);

  // Now, we calculate the x->r transform
  xgrid.Init(0.0, xMax, NumPoints);
  Temp(0) = 0.0;
  IntegrateFirstOrderNS(xgrid, 0, NumPoints-1, Temp, dr_dx, pot);
  StartDeriv = sqrt(A0);
  EndDeriv = 1.0;

  r_of_x.Init(&xgrid, Temp, StartDeriv, EndDeriv);
}
