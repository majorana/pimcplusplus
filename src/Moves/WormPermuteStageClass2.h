#ifndef WORM_PERMUTE_STAGE2_CLASS_H
#define WORM_PERMUTE_STAGE2_CLASS_H

#include "PermuteStageClass.h"


class WormPermuteStage2Class : public PermuteStageClass
{
private:
  ObservableVecDouble1 AcceptanceRatioVar;
  ObservableVecInt1 AcceptanceTotalVar;
  Array<int,1> NumAccepted;
  Array<int,1> NumAttempted;
  int CountPtclInOpenLoop();
  bool PtclInOpenLoop(int checkPtcl);
  bool OnlyX;
  //This is the list of the particle to which you can permute onto.
  list<int> PermuteOntoList;
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
  void Accept ();
  void Reject();
  void WriteRatio();
  WormPermuteStage2Class (PathDataClass &pathData, int speciesNum, int numLevels,IOSectionClass &outSection) : 
    PermuteStageClass(pathData, speciesNum, numLevels,outSection),
    AcceptanceRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator),    AcceptanceTotalVar("Perms Tried",OutSection,PathData.Path.Communicator)

  {
    NumAccepted.resize(4);
    NumAttempted.resize(4);
    NumAccepted=0;
    NumAttempted=0;
  }
};



#endif
