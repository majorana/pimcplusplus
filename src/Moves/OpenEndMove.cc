#include "OpenEndMove.h"
#include "BisectionStageClass.h"
#include "EndStageClass.h"

void OpenEndMoveClass::WriteRatio()
{

  //Do nothing for now
}

void OpenEndMoveClass::Read(IOSectionClass &in)
{
  string speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("name", Name));
  SpeciesNum = PathData.Path.OpenSpeciesNum;
  StageClass* endStage;
  endStage=new EndStageClass(PathData,NumLevels);
  endStage->Read(in);
  Stages.push_back (endStage);
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level,
							     OutSection);
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



