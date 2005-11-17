#ifndef MULTI_STAGE_H
#define MULTI_STAGE_H


#include <list>
#include "MoveBase.h"
#include "../Observables/ObservableBase.h"
#include "../Actions/ActionBase.h"
#include "StageClass.h"

///One method of making a move is to build it out of a set of
///stages. This class allows for such a construct. There is a list of
///stages (in the variable Stages).  By default, this moves iterates
///over the stages, running each one sequentially. Each stage must
///return true or false indicating
// whether or not that stage is accepted. If the stage is
///accepted then the next stage is run. 
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
  virtual void WriteRatio();

  ///Why was this MakeMove()=0 and virtual?
  virtual void MakeMove();
  MultiStageClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    ParticleMoveClass(pathData,outSection) 
  {
    //do nothing for now
  }
};


#endif
