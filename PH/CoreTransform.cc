#include "..//Integration/Integrate.h"
#include "CoreTransform.h"

scalar dx_dr (scalar r, scalar x, void *PHptr)
{
  PseudoHamiltonian &PH = *(PseudoHamiltonian *)PHptr;
  scalar A, B, V, dAdr;
  PH.ABV(r, A, B, V, dAdr);

  return (1.0/sqrt(A));
}

scalar dr_dx (scalar x, scalar r, void *PHptr)
{
  PseudoHamiltonian &PH = *(PseudoHamiltonian *)PHptr;

  scalar A, B, V, dAdr;
  PH.ABV(r, A, B, V, dAdr);

  return (sqrt(A));
}

void
CoreTransform::Initialize(PseudoHamiltonian *PH, int NumPoints)
{

  // First, we calculate the r->x transform by integrating
  rMax = PH->CoreRadius;
  rgrid.Init (0.0, rMax, NumPoints);
  Array<scalar,1> Temp(NumPoints);

  Temp(0) = 0.0;
  IntegrateFirstOrderNS(rgrid, 0, NumPoints-1, Temp, dx_dr, PH);
  xMax = Temp(NumPoints-1);

  scalar A0, B0, V0, dA0;
  PH->ABV(0.0, A0, B0, V0, dA0);
  scalar StartDeriv = 1.0/sqrt(A0);
  scalar EndDeriv = 1.0;

  x_of_r.Init(&rgrid, Temp, StartDeriv, EndDeriv);

  // Now, we calculate the x->r transform
  xgrid.Init(0.0, xMax, NumPoints);
  Temp(0) = 0.0;
  IntegrateFirstOrderNS(xgrid, 0, NumPoints-1, Temp, dr_dx, PH);
  StartDeriv = sqrt(A0);
  EndDeriv = 1.0;

  r_of_x.Init(&xgrid, Temp, StartDeriv, EndDeriv);
}
