#ifndef WORM_SWAP_MOVE_H
#define WORM_SWAP_MOVE_H

#include "../PathDataClass.h"
// #include "WormStage.h"
#include "MultiStage.h"



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class WormSwapMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  //  WormStageClass WormSwapStage;
public:
  ///Specifies how big of a chunk to change while making swap. May
/// eventually change to be more flexible
  int NumLevels;
  void ChooseTimeSlices();
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  //  void MakeMove();
  WormSwapMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection)
    //    WormSwapStage(pathData, outSection)
  {
    DumpFreq = 20;
  }
};


#endif
