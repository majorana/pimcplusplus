#include "PeriodicSpline.h"

void TestPeriodic()
{
  LinearGrid xGrid(0.0, 4.0, 5);
  Array<double,1> vals(5);
  vals(0) = 0.0;
  vals(1) = 1.0;
  vals(2) = 0.0;
  vals(3) = -1.0;
  vals(4) = 0.0;

  PeriodicSpline spline;
  spline.Init (&xGrid, vals);

  for (double x=0.0; x < 4.0; x+= 0.01)
    fprintf (stdout, "%1.16e %1.16e\n", x, spline(x));

}


main()
{
  TestPeriodic();
}
