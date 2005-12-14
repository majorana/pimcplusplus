#ifndef SHORT_RANGE_CLASS_H
#define SHORT_RANGE_CLASS_H

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>

#include "ShortRangeOnClass.h"

/// The ShortRangeClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class ShortRangeClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<bool,1> DoPtcl;
  ShortRangeOnClass ToCheck;
  int TotalTime;
  /// These are the coefficients used for the low-variance estimator
  /// for the gradient
  Array<double,1> ck;
  int NumBasisFuncs, m;
  double Router;
  bool UseLowVariance;
  void Setup_ck();
  inline double g(double r);
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  void GradAction (int slice1, int slice2, const Array<int,1> &ptcls,
		   int level, Array<dVec,1> &gradVec);
  ShortRangeClass (PathDataClass &pathData,
		   Array<PairActionFitClass*, 2> &pairMatrix);
};

#endif
