#include "MultiStage.h"

void StageClass::Accept()
{
  ///Do nothing for now
}

void StageClass::Reject()
{
  ///Do nothing for now
}

void MultiStageClass::Read(IOSectionClass& in)
{

  ///do nothing for now
}
void MultiStageLocalClass::MakeMove()
{
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end() && toAccept){
    SetMode(OLDMODE);
    double oldAction=(*stageIter)->StageAction(Slice1,Slice2,ActiveParticles);
    SetMode(NEWMODE);
    double sampleRatio=(*stageIter)->Sample(Slice1,Slice2,ActiveParticles);    
    double newAction = (*stageIter)->StageAction(Slice1,Slice2,ActiveParticles);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=sampleRatio+currActionChange-prevActionChange;

    
    if (-logAcceptProb<log(PathData.Path.Random.Local())) ///reject conditin
      toAccept=false;
    //    if (toAccept)
    //      
    //      NumAccepted(levelCounter)++;
    //    else
    //      NumRejected(levelCounter)++;
    prevActionChange=currActionChange;
    stageIter++;
  }
  if (toAccept)
    Accept();
  else 
    Reject();
}

void MultiStageClass::Accept()
{
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();
       stageIter!=Stages.end();stageIter++){
    (*stageIter)->Accept();
  }  
}

void MultiStageClass::Reject()
{
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator 
	 stageIter=Stages.begin();stageIter!=Stages.end();stageIter++){
    (*stageIter)->Reject();
  }
  
}
void MultiStageCommonClass::MakeMove()
{
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end() && toAccept){
    SetMode(OLDMODE);
    double oldAction=(*stageIter)->StageAction(Slice1,Slice2,ActiveParticles);
    SetMode(NEWMODE);
    double sampleRatio=(*stageIter)->Sample(Slice1,Slice2,ActiveParticles);    
    double newAction = (*stageIter)->StageAction(Slice1,Slice2,ActiveParticles);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=sampleRatio+currActionChange-prevActionChange;

    
    if (-logAcceptProb<log(PathData.Path.Random.Common())) ///reject conditin
      toAccept=false;
    //    if (toAccept)
    //      
    //      NumAccepted(levelCounter)++;
    //    else
    //      NumRejected(levelCounter)++;
    prevActionChange=currActionChange;
    stageIter++;
  }
  if (toAccept){
    PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  }
  else {
    PathData.RejectMove(Slice1,Slice2,ActiveParticles);
  }
}
