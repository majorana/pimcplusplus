#include "ColorMap.h"

void
ColorMap::Init (double min, double max, ColorMapType map)
{
  Min = min;
  Max = max;

  Array<double,1> r, g, b, a;
  if (map == BLUE_WHITE_RED) {
    r.resize(10); g.resize(10); b.resize(10); a.resize(10);
    r = 0.00, 0.00, 0.00, 0.10, 0.50, 0.80, 0.85, 0.85, 0.70, 0.50;
    g = 0.00, 0.05, 0.20, 0.50, 0.85, 0.85, 0.50, 0.20, 0.10, 0.05;
    b = 0.50, 0.80, 0.85, 0.85, 0.80, 0.50, 0.20, 0.05, 0.00, 0.00;
    a = 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99;
  }
  else {
    cerr << "Unknown map type in ColorMap::Init.\n";
    abort();
  }
  BoundaryCondition<double> fBC(FLAT), nBC(NATURAL);
  Splines[0].Init (min, max, r, true, fBC, nBC);
  Splines[1].Init (min, max, g, true, fBC, fBC);
  Splines[2].Init (min, max, b, true, nBC, fBC);
  Splines[3].Init (min, max, a, true, nBC, nBC);
//   FILE *fout = fopen ("colormap.dat", "w");
//   for (double x=min; x<=max; x+=0.001)
//     fprintf (fout, "%10.6e %10.6e %10.6e %10.6e %10.6e\n",
// 	     x, Splines[0](x), Splines[1](x), Splines[2](x), Splines[3](x));
//   fclose (fout);
}
