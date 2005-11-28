#include "ShortRangePotClass.h"
#include "../PathDataClass.h"

ShortRangePotClass::ShortRangePotClass 
(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix) :
  PotentialBaseClass (pathData), PairMatrix(pairMatrix)
{
  // Do nothing 
}

double 
ShortRangePotClass::V(int slice)
{
  double val = 0.0;
  for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      double dist;
      dVec disp;
      
      PathData.Path.DistDisp(slice, ptcl1, ptcl2, dist, disp);
      PairActionFitClass& pa=
	*(PairMatrix(species1, PathData.Path.ParticleSpeciesNum(ptcl2)));
      
      val += pa.V(dist);
      if (pa.IsLongRange() && PathData.Actions.UseLongRange)
	val -= pa.Vlong(dist);
    }
  }
  return val;
}
