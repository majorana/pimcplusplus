#include "PermuteStageClass.h"

double NoPermuteStageClass::Sample (int &slice1, int &slice2,
				    Array<int,1> &activeParticles)
{
  if (activeParticles(0) == -1)
    activeParticles.resize (1);
  activeParticles(0) = ChooseParticle();
  return 1.0;
}

int NoPermuteStageClass::ChooseParticle()
{
  return (PathData.Path.Random.LocalInt 
    (PathData.Path.Species(SpeciesNum).NumParticles));
}
