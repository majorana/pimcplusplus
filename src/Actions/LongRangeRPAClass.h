#ifndef LONG_RANGE_RPA_CLASS_H
#define LONG_RANGE_RPA_CLASS_H

#include "ActionBase.h"
#include "../Common/PairAction/PAFit.h"

/// LongRangeRPAClass is the Random Phase Approximation corrected
/// version of LongRangeClass.  This class solves the couple
/// differential equations derived from the random phase approximation
/// to improve the converge of the actions as \f$ \tau \rightarrow 0 \f$.
class LongRangeRPAClass : public ActionBaseClass
{
protected:
  Array<PairActionFitClass*,2> &PairMatrix;
  Array<PairActionFitClass*,1> &PairArray;
  LinearGrid LongGrid;

  /// This must be called after all of the OptimizedBreakup_x's
  Array<double,1> Integrand(double t, const Array<double,1> &Uvec);

  void Test();
  bool TaskIsU;
  /// These keep track of which level and k vector we are working on
  /// when solving the RPA equations.
  int Level, ki;
public:
  inline Array<double,1> operator()(double t, Array<double,1> uwvec)
  { return Integrand(t, uwvec); }

public:
  bool UseBackground;
  /// NOTE: this can only be called after LongRangeClass's Init has
  /// been called.
  void Init(IOSectionClass &in);
  double Action (int slice1, int slice2, 
		   const Array<int,1> &activeParticles, int level);
  double d_dBeta(int slice1, int slice2, int level);
  LongRangeRPAClass (PathDataClass &pathData,
		     Array<PairActionFitClass*, 2> &pairMatrix,
		     Array<PairActionFitClass*, 1> &pairArray);

};

#endif
