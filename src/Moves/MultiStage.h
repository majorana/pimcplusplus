#ifndef MULTI_STAGE_H
#define MULTI_STAGE_H


#include "list.h"
#include "MoveClass.h"
#include "../Actions/ActionBase.h"

class StageClass
{
protected:
  PathDataClass &PathData;
public:
  int BisectionLevel;
  list<ActionBase*> Actions;
  ///The highest stage will set the slices and activeParticles
  ///This returns transition probability T(new->old)/T(old->new)
  virtual double Sample (int &slice1,int &slice2,
			 Array<int,1> activeParticles); 
  double StageAction(int startSlice,int endSlice,
		     const Array<int,1> &changedParticles);
  StageClass(PathDataClass &pathData) :PathData(pathData)
  {
    //do nothing for now
  }
};


double StageClass::StageAction(int startSlice,int endSlice,
		   const Array<int,1> &changedParticles)
{
  double TotalAction=0.0;
  list<ActionBase*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()){
    TotalAction += 
      ((*actionIter)->Evaluate(startSlice, endSlice, changedParticles,
			       BisectionLevel);
    actionIter++;
  }
  return TotalAction;
}


class MultiStageClass : public ParticleMoveClass
{
protected:
  list<StageClass*> Stages;
  int NumSteps;
public:
  void Read(IOSectionClass &io);
  void MakeMove();
  MultiStageClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    ParticleMoveClass(pathData,outSection) 
  {
    //do nothing for now
  }
};

void MultiStageClass::MakeMove()
{

  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end()){
    SetMode(OLDMODE);
    double oldAction=stageIter->StageAction(slice1,slice2,activeParticles);
    SetMode(NEWMODE);
    double sampleRatio=stageIter->Sample(slice1,slice2,activeParticles);    
    double newAction = stageIter->StageAction(slice1,slice2,activeParticles)
    double currActionChange=newAction-oldAction;
    double logAcceptProb=sampleRatio+currActionChange-prevActionChange;
    if (-logAcceptProb<log(PathData.Path.Random.Local())) ///reject conditin
      toAccept=false;
    if (toAccept)
      NumAccepted(levelCounter)++;
    else
      NumRejected(levelCounter)++;
    prevActionChange=currActionChange;
    stageIter++;
  }
}


#endif
