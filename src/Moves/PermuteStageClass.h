#ifndef PERMUTE_STAGE_CLASS_H
#define PERMUTE_STAGE_CLASS_H

#include "MultiStage.h"

class PermuteStageClass : public LocalStageClass
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
  virtual bool Attempt (int &slice1, int &slice2, 
			Array<int,1> &activeParticles, double &prevActionChange) = 0;
  virtual void InitBlock();
  virtual void Read (IOSectionClass &in);
  virtual void Accept();
  virtual void Reject();
  PermuteStageClass(PathDataClass &pathData, int speciesNum,
		    int numLevels) : 
    LocalStageClass (pathData), 
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
  bool Attempt (int &slice1, int &slice2, 
		 Array<int,1> &activeParticles, double &prevActionChange);
  NoPermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels) 
    : PermuteStageClass(pathData, speciesNum, numLevels)
  {
    // do nothing for now
  }
};


class TablePermuteStageClass : public PermuteStageClass
{
private:
  //  PermuteTableClass Table1, Table2;
  //  PermuteTableClass *Forw, *Rev;
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
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);
  bool Attempt (int &slice1, int &slice2, 
		   Array<int,1> &activeParticles, double &prevActionChange);
  TablePermuteStageClass (PathDataClass &pathData, 
			  int speciesNum, int numLevels) : 
    PermuteStageClass(pathData, speciesNum, numLevels)
  {
    // do nothing for now
  }
};


#endif
