#ifndef SHORT_RANGE_POT_CLASS
#define SHORT_RANGE_POT_CLASS

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>

class ShortRangePotClass : public PotentialBaseClass
{
private:
  Array<PairActionFitClass*,2> &PairMatrix;
public:
  double V (int slice);

  ShortRangePotClass (PathDataClass &pathData,
		      Array<PairActionFitClass*,2> &pairMatrix);
};



#endif
