#include "OpenEndMove.h"
#include "BisectionStage.h"
#include "EndStage.h"

void OpenEndMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();
  //Do nothing for now
}

void OpenEndMoveClass::MakeMove()
{
  //  cerr<<"GRR!!  I'm HERE!"<<endl;
  //  sleep(1000);
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }
  if (toAccept)
    MultiStageClass::Accept();
  else 
    MultiStageClass::Reject();
  //GRR!  MultiStageClass::MoveClass::MakeMove();
}

void OpenEndMoveClass::Read(IOSectionClass &in)
{
  string speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("name", Name));
  SpeciesNum = PathData.Path.OpenSpeciesNum;
  StageClass* endStage;
  endStage=new EndStageClass(PathData,NumLevels,IOSection);
  endStage->Read(in);
  Stages.push_back (endStage);
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = 
      new BisectionStageClass (PathData, level, IOSection);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    if (level==0)
      newStage->Actions.push_back(&PathData.Actions.OpenLoopImportance);
    if (level>0)
      newStage->Actions.push_back(&PathData.Actions.ShortRangeApproximate);
    else
      newStage->Actions.push_back(&PathData.Actions.ShortRange);
    if (level == 0) {
      if (PathData.Path.LongRange){
	if (PathData.Actions.UseRPA)
	  newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
      ///If it's David's long range class then do this
	else if (PathData.Path.DavidLongRange){
	  newStage->Actions.push_back(&PathData.Actions.DavidLongRange);
	}
	////
	else
	  newStage->Actions.push_back(&PathData.Actions.LongRange);

	
      }
      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	cerr << "Adding fermion node action for species " 
	     << speciesName << endl;

	newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
    }

    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }

}



