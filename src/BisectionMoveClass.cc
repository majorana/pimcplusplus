#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "Common.h"


void ShiftMove::makeMove()
{//Remember to mark Actions dirty!!!
  int numTimeSlicesToShift=(int)floor(sprng()*PathData->NumTimeSlices);
  for (int counter=0;counter<PathData->IdenticalParticleArray.size();counter++)
    {
      (PathData->IdenticalParticleArray)(counter).Path.shiftData(numTimeSlicesToShift,PathData->Communicator);
    }
}

    

BisectionMoveClass::BisectionMoveClass()
{

}

void BisectionMoveClass::makeMove()
{
  bool toAccept=true;
  double oldLogSampleProb;
  double newLogSampleProb;
  Array<int,1> theParticles;
  double logSampleProb;
  int EndTimeSlice=1<<NumLevels+StartTimeSlice;
  double prevActionChange=0;


  ChooseParticles();   
  for (int levelCounter=NumLevels;levelCounter>0;levelCounter--){

    setMode(OLDMODE);
    bool toAccept=true;
    double oldAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,StartTimeSlice,EndTimeSlice,levelCounter-1);
    setMode(NEWMODE);
    (*PathData).TotalAction.SampleParticles(ActiveParticles,StartTimeSlice,EndTimeSlice,levelCounter,newLogSampleProb,oldLogSampleProb);
    double newAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,StartTimeSlice,EndTimeSlice, levelCounter-1);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=oldLogSampleProb-newLogSampleProb+currActionChange-prevActionChange;
    if (-logAcceptProb<log(sprng())){///reject conditin
      toAccept=false;
      break;

    }

    prevActionChange=currActionChange;
    
  }
  if (toAccept ==true ){
    (*PathData).acceptMove(ActiveParticles,StartTimeSlice,EndTimeSlice);
  }
  else {
    (*PathData).rejectMove(ActiveParticles,StartTimeSlice,EndTimeSlice);
  }
    

}
