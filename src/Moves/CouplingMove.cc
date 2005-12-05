#include "CouplingMove.h"
#include "BisectionStage.h"
#include "CouplingStage.h"

void CouplingMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();
  //Do nothing for now
}

void CouplingMoveClass::MakeMove()
{

  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  int stageCount=0;
  while (stageIter!=Stages.end() && toAccept){
    //    cerr<<"This is stage "<<stageCount<<endl;
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
    //    cerr<<endl;
  }
  if (toAccept)
    MultiStageClass::Accept();
  else 
    MultiStageClass::Reject();
}

void CouplingMoveClass::Read(IOSectionClass &in)
{
  StageClass* coupleStage;
  coupleStage=new CouplingStageClass(PathData,NumLevels,IOSection);
  coupleStage->Read(in);
  coupleStage->Actions.push_back(&PathData.Actions.Kinetic);
  coupleStage->Actions.push_back(&PathData.Actions.ShortRange);
  Stages.push_back (coupleStage);
}



