#include "SplinePot.h"

double SplinePot::V(double r)
{
  return Spline(r);
}

double SplinePot::dVdr(double r)
{
  return Spline.Deriv(r);
}

double SplinePot::d2Vdr2(double r)
{
  return Spline.Deriv2(r);
}

void SplinePot::Read(IOSectionClass &in)
{
  assert(in.OpenSection("Grid"));
  Grid *grid = ReadGrid(in);
  Array<double,1> data;
  assert(in.ReadVar("SplineData", data));
  Spline.Init (grid, data);
}


void SplinePot::Write(IOSectionClass &out)
{
  out.WriteVar ("Type", "Spline");
  out.NewSection("Grid");
  Spline.grid->Write(out);
  out.CloseSection();

  out.WriteVar ("SplineData", Spline.Data());
}
