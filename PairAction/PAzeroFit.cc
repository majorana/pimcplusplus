#include "PAFit.h"

double PAzeroFitClass::U(double q, double z, double s2, int level)
{
  return 0.0;
}

double PAzeroFitClass::dU(double q, double z, double s2, int level)
{
  return 0.0;
}

double PAzeroFitClass::V(double r)
{
  return 0.0;
}

bool PAzeroFitClass::Read (IOSectionClass &in,
				double smallestBeta, int NumBetas)
{
  //  SmallestBeta = smallestBeta;
  // Read Particles;
    assert(in.OpenSection("Fits"));
    assert(in.OpenSection("Particle1"));
    Particle1.Read(in);
    in.CloseSection();
    assert(in.OpenSection("Particle2"));
    Particle2.Read(in);
    in.CloseSection();
  //  lambda = Particle1.lambda + Particle2.lambda;
  //  assert (lambda == 0.0);

  // Read Potential;
  //  assert(in.OpenSection("Potential"));
  //  Potential = ReadPH(in);
  //  in.CloseSection();
  //  in.CloseSection();
  return true;
}



bool PAzeroFitClass::IsLongRange()
{
  return false;
}
