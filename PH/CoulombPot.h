#ifndef COULOMB_POT_H
#define COULOMB_POT_H

#include "PotentialBase.h"

class CoulombPot : public Potential
{
public:
  double Z1Z2;
  
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);
};

#endif
