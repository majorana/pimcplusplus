#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "CoulombPot.h"
#include "ScreenedPot.h"
#include "QuinticPH.h"

inline Potential* ReadPotential (IOSectionClass &in)
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
  else {
    cerr << "Unrecognize potential type \"" << type << "\".  Exitting.\n";
    exit(1);
  }
  pot->Read(in);
  return pot;
}

#endif
