#include "ToppHopfield.h"

double 
ToppHopfieldPot::V(double r)
{
  if (r < rc)
    return (V0*cos(a*r)+b);
  else
    return (-Z/r);
}

double
ToppHopfieldPot::dVdr(double r)
{
  if (r < rc)
    return (-V0*a*sin(a*r));
  else
    return (Z/(r*r));
}

double
ToppHopfieldPot::d2Vdr2(double r)
{
  if (r < rc)
    return (-V0*a*a*cos(a*r));
  else
    return (-2.0*Z/(r*r*r));
}

void
ToppHopfieldPot::Read(IOSectionClass &in)
{
  assert (in.ReadVar("Z",  Z));
  assert (in.ReadVar("V0", V0));
  assert (in.ReadVar("a",  a));
  assert (in.ReadVar("b",  b));
  assert (in.ReadVar("rc", rc));
}

void
ToppHopfieldPot::Write(IOSectionClass &out)
{
  out.WriteVar("Type", "ToppHopfield");
  out.WriteVar("Z",  Z);
  out.WriteVar("V0", V0);
  out.WriteVar("a",  a);
  out.WriteVar("b",  b);
  out.WriteVar("rc", rc);
}
