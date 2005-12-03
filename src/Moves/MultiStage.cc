#include "MultiStage.h"


void MultiStageClass::Read(IOSectionClass& in)
{
  ///do nothing for now
}

void MultiStageClass::WriteRatio()
{
   list<StageClass*>::iterator stageIter=Stages.begin();
   double prevActionChange=0.0;
   while (stageIter!=Stages.end()){
     //    cerr<<"Some stage is writing their ratio"<<endl;
     (*stageIter)->WriteRatio();
     stageIter++;
   }  
   MoveClass::WriteRatio();
}


void MultiStageClass::MakeMove()
{
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;

  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }

  if (toAccept)
    Accept();
  else 
    Reject();
  MoveClass::MakeMove();
}


void MultiStageClass::Accept()
{
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();
       stageIter!=Stages.end();stageIter++){
    (*stageIter)->Accept();
  }  
  NumAccepted++;
}

void MultiStageClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator 
	 stageIter=Stages.begin();stageIter!=Stages.end();stageIter++){
    (*stageIter)->Reject();
  }
}

