#include "SplinePot.h"

double SplinePot::V(double r)
{
  if (r <= Spline.grid->End)
    return Spline(r);
  else {
#ifdef BZ_DEBUG
    cerr << "r outside grid in SplinePot:  " << r << endl;
#endif
    return 0.0;
  }
}

double SplinePot::dVdr(double r)
{
  if (r <= Spline.grid->End)
    return Spline.Deriv(r);
  else {
#ifdef BZ_DEBUG
    cerr << "r outside grid in SplinePot:  " << r << endl;
#endif
    return 0.0;
  }
}

double SplinePot::d2Vdr2(double r)
{
  if (r <= Spline.grid->End)
    return Spline.Deriv2(r);
  else {
#ifdef BZ_DEBUG
    cerr << "r outside grid in SplinePot:  " << r << endl;
#endif
    return 0.0;
  }
}

void SplinePot::Read(IOSectionClass &in)
{
  assert(in.OpenSection("Grid"));
  Grid *grid = ReadGrid(in);
  in.CloseSection(); // "Grid" 
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
