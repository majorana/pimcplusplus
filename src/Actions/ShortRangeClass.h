#ifndef SHORT_RANGE_CLASS_H
#define SHORT_RANGE_CLASS_H

#include "ActionBase.h"
#include "../Common/PairAction/PAFit.h"

class ShortRangeClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<bool,1> DoPtcl;
public:
  void Read (IOSectionClass &in);
  double Evaluate (int slice1, int slice2, 
		   const Array<int,1> &activeParticles, int level);
  ShortRangeClass (PathDataClass &pathData,
		   Array<PairActionFitClass*, 2> &pairMatrix);
};

#endif
