#ifndef SHORT_RANGE_ON_CLASS_H
#define SHORT_RANGE_ON_CLASS_H

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>

/// The ShortRangeClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class ShortRangeOnClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<bool,1> DoPtcl;
public:
  void Read (IOSectionClass &in);
  double Action (int slice1, int slice2, 
		 const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  ShortRangeOnClass (PathDataClass &pathData,
		   Array<PairActionFitClass*, 2> &pairMatrix);
};

#endif
