#include "CoulombPot.h"

double CoulombPot::V (double r)
{
  return (Z1Z2 / r);
}

double CoulombPot::dVdr (double r)
{
  return (-Z1Z2/(r*r));
}

double CoulombPot::d2Vdr2 (double r)
{
  return (2.0*Z1Z2/(r*r*r));
}

void CoulombPot::Write(IOSectionClass &out)
{
  out.WriteVar("Z1Z2", Z1Z2);
}

void CoulombPot::Read(IOSectionClass &in)
{
  assert(in.ReadVar("Z1Z2", Z1Z2));
}
