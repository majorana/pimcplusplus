#include "Potential.h"

Potential* ReadPotential (IOSectionClass &in)
{
  string type;
  Potential *pot;
  assert (in.ReadVar ("Type", type));
  if (type == "Coulomb")
    pot = new CoulombPot;
  else if (type == "QuinticPH")
    pot = new QuinticPH;
  else if (type == "Screened")
    pot = new ScreenedPot;
  else if (type == "Spline")
    pot = new SplinePot;
  else if (type == "HeAziz")
    pot = new HeAzizPot;
  else {
    cerr << "Unrecognize potential type \"" << type << "\".  Exitting.\n";
    exit(1);
  }
  pot->Read(in);
  return pot;
}
