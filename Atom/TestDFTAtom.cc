#include "NewAtom.h"

void TestDFTAtom()
{
  CoulombPot barePot;
  barePot.Z1Z2 = -1.0;
  OptimalGrid grid(1.0, 50.0);
  DFTAtom atom;
  atom.RadialWFs.resize(1);
  atom.RadialWFs(0).l = 0;
  atom.RadialWFs(0).Occupancy = 1.0;
  atom.RadialWFs(0).Energy = -0.3;
  atom.SetGrid (&grid);
  atom.SetBarePot (&barePot);

  atom.newMix = 0.9;
  atom.SolveSC();
}



main()
{
  TestDFTAtom();
}
  
  
