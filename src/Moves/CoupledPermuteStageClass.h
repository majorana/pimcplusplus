#ifndef COUPLED_PERMUTE_STAGE_CLASS_H
#define COUPLED_PERMUTE_STAGE_CLASS_H

#include "MultiStage.h"
#include "PermuteStageClass.h"
#include "PermuteTableClass.h"


class CoupledPermuteStageClass : public PermuteStageClass
{
private:
  PermuteTableClass Table1, Table2;
  PermuteTableClass *Forw, *Rev;
public:
  /// This function will construct a new permutation if
  /// activeParticles is set to the array, [ -1 ];  In this case,
  /// it will set the activeParticles.  It returns 1.0
  /// as the transition probability ratio at this time.
  /// If called with a valid set of particles in activeParticles,
  /// it changes nothing, but returns the transition probability
  /// ratio for the sampling.  This is so we can avoid calculating
  /// that ratio if the move is rejected, saving time.  Thus, this
  /// function is called twice during a successful multistage move.

  void InitBlock();
  void Read (IOSectionClass &in);
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);
  bool Attempt (int &slice1, int &slice2, 
		   Array<int,1> &activeParticles, double &prevActionChange);
  CoupledPermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels) : 
    PermuteStageClass(pathData, speciesNum, numLevels),
    Table1(pathData), Table2(pathData)
  {
    Forw = &Table1;
    Rev  = &Table2;
  }
};


#endif
