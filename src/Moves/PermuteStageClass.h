#ifndef PERMUTE_STAGE_CLASS_H
#define PERMUTE_STAGE_CLASS_H

#include "MultiStage.h"

class PermuteStageClass : public StageClass
{
protected:
  int SpeciesNum, NumLevels;
public:
  /// This function is responsible for setting the activeParticles if
  /// activeParticles(0) == -1; It is called a second time after the
  /// bisection stages.  We impose the condtion that the product of
  /// the two return values must be equal to true transition
  /// probablity ratio.  This allows avoiding computing the reverse
  /// probability if the move is rejected before this stage.
  virtual double Sample (int &slice1, int &slice2, 
			 Array<int,1> &activeParticles) = 0;
  virtual void InitBlock();
  virtual void Read (IOSectionClass &in);
  virtual void Accept();
  virtual void Reject();
  PermuteStageClass(PathDataClass &pathData, int speciesNum,
		    int numLevels) : 
    StageClass (pathData), 
    SpeciesNum (speciesNum), 
    NumLevels(numLevels)
  {
    // do nothing for now 
  }
};

class NoPermuteStageClass : public PermuteStageClass
{
private:
  int ChooseParticle();
  int NumLevels;
  int NumMoves, NumAccepted;
public:
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);
  NoPermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels) 
    : PermuteStageClass(pathData, speciesNum, numLevels)
  {
    // do nothing for now
  }
};


class TablePermuteStageClass : public PermuteStageClass
{
public:
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);

  TablePermuteStageClass (PathDataClass &pathData, 
			  int speciesNum, int numLevels) : 
    PermuteStageClass(pathData, speciesNum, numLevels)
  {
    // do nothing for now
  }
};


#endif
