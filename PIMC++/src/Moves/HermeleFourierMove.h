#ifndef HERMELE_FOURIER_MOVE_H
#define HERMELE_FOURIER_MOVE_H

#include "../PathDataClass.h"
#include "MultiStage.h"


/// This stage attempts to displace a list of whole paths.  It should
/// only be used for non-permuting particles.
class HermeleFourierStageClass : public CommonStageClass
{
public:
  /// This is the width of the gaussian distribution for the
  /// displacement vector.
  double Sigma;
  /// This does the actual displacement of the path.  All processors
  /// within a single close must displace by the same amount.
  double Sample (int &slice1, int &slice2,
		 Array <int,1> &activeParticles);
  HermeleFourierStageClass (PathDataClass &pathData,IOSectionClass &outSection) :
    CommonStageClass (pathData,outSection) 
  {
    // Do nothing for now.
  }
};



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class HermeleFourierMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  HermeleFourierStageClass HermeleFourierStage;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  HermeleFourierMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    HermeleFourierStage(pathData, outSection)
  {
    DumpFreq = 20;
  }
};


#endif
