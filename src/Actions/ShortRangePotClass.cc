#include "ShortRangePotClass.h"

ShortRangePotClass::ShortRangePotClass 
(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix) :
  PotentialBaseClass (pathData), PairMatrix(pairMatrix)
{
  // Do nothing 
}

double ShortRangePotClass::V(int slice)
{
  return 0.0;
}
