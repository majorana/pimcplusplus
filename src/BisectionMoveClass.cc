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
  int EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
  double prevActionChange=0;
  
  //  cerr<<"At the beginning fo the makeMove the size is ";
  //  cerr <<  PathData->IdenticalParticleArray.size()<<endl;
  ChooseParticles();   
  //////////  for (int levelCounter=NumLevels;levelCounter>0;levelCounter--){
  int levelCounter=NumLevels-1;

  while (levelCounter>=0 && toAccept==true){
    //    cerr<<"At the level Counter in Bisection makeMove being  "<<levelCounter<<" the size is ";
    //    cerr <<  PathData->IdenticalParticleArray.size()<<endl;
    setMode(OLDMODE);
    toAccept=true;
    double oldAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,StartTimeSlice,EndTimeSlice,levelCounter);
    oldLogSampleProb = (*PathData).TotalAction.LogSampleProb(ActiveParticles,StartTimeSlice,EndTimeSlice,levelCounter);
    setMode(NEWMODE);
    newLogSampleProb = (*PathData).TotalAction.SampleParticles(ActiveParticles,StartTimeSlice,EndTimeSlice,levelCounter);
    //cerr << "newLogSampleProb = " << newLogSampleProb << endl;
    //cerr << "oldLogSampleProb = " << oldLogSampleProb << endl;
    double newAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,StartTimeSlice,EndTimeSlice, levelCounter);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=oldLogSampleProb-newLogSampleProb-currActionChange+prevActionChange;
    //cerr << "prevActionChange = " << prevActionChange << endl;
    //cerr << "logAcceptProb = " << logAcceptProb << endl;
    //cerr<<"My new action is "<<newAction<<" and my old action was "<<oldAction<<endl;
    if (-logAcceptProb<log(sprng())){///reject conditin
      toAccept=false;
      //      break;

    }

    prevActionChange=currActionChange;
    levelCounter--;
  }
  if (toAccept ==true ){
    (*PathData).acceptMove(ActiveParticles,StartTimeSlice,EndTimeSlice);
    //cout<<"I'm accepting! I'm accepting!"<<endl;
  }
  else {
    (*PathData).rejectMove(ActiveParticles,StartTimeSlice,EndTimeSlice);
    //  cout<<"I'm rejecting! I'm rejecting!"<<endl;
  }
    

}
