#ifndef LONG_RANGE_POT_CLASS
#define LONG_RANGE_POT_CLASS

#include "ActionBase.h"
#include "../Common/PairAction/PAFit.h"

class LongRangePotClass : public PotentialBaseClass
{
private:
  Array<PairActionFitClass*,2> &PairMatrix;
public:
  double V (int slice);

  LongRangePotClass (PathDataClass &pathData,
		     Array<PairActionFitClass*,2> &pairMatrix);
};



#endif
