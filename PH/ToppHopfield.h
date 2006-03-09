#ifndef TOPP_HOPFIELD_H
#define TOPP_HOPFIELD_H

#include "PotentialBase.h"

/// References Phys. Rev. B 55:23 15515 (1997)
///            Phys. Rev. B 7, 1295 (1973)

class ToppHopfieldPot : public Potential
{
public:
  double V0, a, b, rc, Z;

  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);

  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);
};

#endif
