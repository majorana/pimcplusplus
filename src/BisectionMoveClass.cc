#include "BisectionMoveClass.h"
#include "Common.h"

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

    setMode(OLD);
    bool toAccept=true;
    double oldAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,levelCounter-1);
    setMode(NEW);
    (*PathData).TotalAction.SampleParticles(ActiveParticles,StartTimeSlice,EndTimeSlice,levelCounter,newLogSampleProb,oldLogSampleProb);
    double newAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,levelCounter-1);
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
