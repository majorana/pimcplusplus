#include "RadialWF.h"

void TestRadialWF1()
{
  CoulombPot pot;
  pot.Z1Z2 = -1.0;

  OptimalGrid grid(1.0, 50.0);
  RadialWF wf;
  wf.l = 1;
  wf.TotalNodes = 0;
  wf.Energy = -0.25;
  wf.Occupancy = 1.0;
  wf.SetPotential (&pot);
  wf.SetGrid (&grid);
  
  wf.Solve (1.0e-14);
  fprintf (stderr, "Energy = %1.12f\n", wf.Energy);
}



main()
{
  TestRadialWF1();
}
