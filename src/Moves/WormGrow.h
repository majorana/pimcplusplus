#ifndef WORM_GROW_MOVE_H
#define WORM_GROW_MOVE_H

#include "WormStage.h"
#include "MultiStage.h"
#include "../PathDataClass.h"



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class WormGrowMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  WormStageClass WormGrowStage;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  int MaxGrowth;
  WormGrowMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    WormGrowStage(pathData, outSection)
  {
    DumpFreq = 20;
  }
};


#endif
