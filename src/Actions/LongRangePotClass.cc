#include "LongRangePotClass.h"

LongRangePotClass::LongRangePotClass 
(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix) :
  PotentialBaseClass (pathData), PairMatrix(pairMatrix)
{
  // Do nothing 
}

double LongRangePotClass::V(int slice)
{
  return 0.0;
}
