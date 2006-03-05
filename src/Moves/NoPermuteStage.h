#ifndef NO_PERMUTE_STAGE_CLASS_H
#define NO_PERMUTE_STAGE_CLASS_H

#include "PermuteStage.h"
#include "MultiStage.h"
#include "PermuteTableClass.h"
#include "../Observables/ObservableVar.h"

class NoPermuteStageClass : public PermuteStageClass
{
private:
  int ChooseParticle();
  int NumLevels;
  int NumMoves, NumAccepted;
public:
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);
  bool Attempt (int &slice1, int &slice2, 
		 Array<int,1> &activeParticles, double &prevActionChange);

  void NoPermuteStageClass::InitBlock(int &slice1,int &slice2);
  NoPermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels,
		       IOSectionClass &outSection) 
    : PermuteStageClass(pathData, speciesNum, numLevels,outSection)
  {
    // do nothing for now
  }
};




#endif
