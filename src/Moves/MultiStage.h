
#ifndef MULTI_STAGE_H
#define MULTI_STAGE_H


#include "list.h"
#include "MoveBase.h"
#include "../Actions/ActionBase.h"

class StageClass
{
protected:
  PathDataClass &PathData;
public:
  int NumAccepted, NumAttempted;
  int BisectionLevel;
  list<ActionBaseClass*> Actions;
  ///The highest stage will set the slices and activeParticles
  ///This returns transition probability T(new->old)/T(old->new)
  virtual double Sample (int &slice1,int &slice2,
			 Array<int,1> &activeParticles) = 0; 
  virtual void Read (IOSectionClass &in);
  virtual void Accept();
  virtual void Reject();
  inline double StageAction(int startSlice,int endSlice,
			    const Array<int,1> &changedParticles);
  inline double AcceptRatio () 
  { return (double)NumAccepted / (double) NumAttempted; }

  StageClass(PathDataClass &pathData) :
    PathData(pathData), NumAccepted(0), NumAttempted(0)
  {
    // Do nothing for now
  }
};


inline double StageClass::StageAction(int startSlice,int endSlice,
		   const Array<int,1> &changedParticles)
{
  double TotalAction=0.0;
  list<ActionBaseClass*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()){
    TotalAction += 
      ((*actionIter)->Action(startSlice, endSlice, changedParticles,
			     BisectionLevel));
    actionIter++;
  }
  return TotalAction;
}


class MultiStageClass : public ParticleMoveClass
{
protected:
  list<StageClass*> Stages;
  int NumSteps;
  int Slice1,Slice2;
public:
  void Read(IOSectionClass &io);
  void Accept();
  void Reject();

  virtual void MakeMove()=0;
  MultiStageClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    ParticleMoveClass(pathData,outSection) 
  {
    //do nothing for now
  }
};

class MultiStageLocalClass : public MultiStageClass
{
public:
  void MakeMove();
  MultiStageLocalClass(PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData,outSection)
  {
    //do nothing for now
  }
};


class MultiStageCommonClass : public MultiStageClass
{
public:
  void MakeMove();
  MultiStageCommonClass(PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData,outSection)
  {
    //do nothing for now
  }
};


#endif
