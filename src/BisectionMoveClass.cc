#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "Common.h"
#include "SpeciesClass.h"




void ShiftMoveClass::MakeMove()
{//Remember to mark Actions dirty!!!
  //int numTimeSlicesToShift=(int)floor(sprng()*PathData->NumTimeSlices);
  int numTimeSlicesToShift = 5;
  
  if (numTimeSlicesToShift > 0){
    for (int counter=0;counter<PathData.NumSpecies();counter++){
      //      PathData.MoveJoin(1);
      PathData.ShiftData(numTimeSlicesToShift);
      //      PathData.Join=1+numTimeSlicesToShift;
    }
  }
  else if (numTimeSlicesToShift<=0){
    for (int counter=0;counter<PathData.NumSpecies();counter++){
      //      PathData.MoveJoin(PathData.NumTimeSlices()-2); //< -1 so you don't overflow and -1 again because you 
                                                              //don't actually want to be at the entire end
      PathData.ShiftData(numTimeSlicesToShift);
      //      PathData.Join=PathData.NumTimeSlices()-2+numTimeSlicesToShift;
    }
  }
  
    
}

    

void BisectionMoveClass::MakeMove()
{
  bool toAccept=true;
  double oldLogSampleProb;
  double newLogSampleProb;
  Array<int,1> theParticles;
  int EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
  double prevActionChange=0;
  
  //  cerr<<"At the beginning fo the makeMove the size is ";
  //  cerr <<  PathData->SpeciesArray.size()<<endl;
  ChooseParticles();   
  //////////  for (int levelCounter=NumLevels;levelCounter>0;levelCounter--){
  int levelCounter=NumLevels-1;

  while (levelCounter>=0 && toAccept==true){
    // cerr<<"At the level Counter in Bisection makeMove being  "
    // <<levelCounter<<" the size is ";
    //    cerr <<  PathData->SpeciesArray.size()<<endl;
    SetMode(OLDMODE);
    toAccept=true;
    double oldAction = PathData.Action.calcTotalAction
      (StartTimeSlice,EndTimeSlice,ActiveParticles,levelCounter);
    oldLogSampleProb = PathData.Action.LogSampleProb
      (StartTimeSlice,EndTimeSlice,ActiveParticles,levelCounter);
    SetMode(NEWMODE);
    newLogSampleProb = PathData.Action.SampleParticles
      (StartTimeSlice,EndTimeSlice,ActiveParticles,levelCounter);
    //cerr << "newLogSampleProb = " << newLogSampleProb << endl;
    //cerr << "oldLogSampleProb = " << oldLogSampleProb << endl;
    double newAction = PathData.Action.calcTotalAction
      (StartTimeSlice,EndTimeSlice, ActiveParticles,levelCounter);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=
      -oldLogSampleProb+newLogSampleProb+currActionChange-prevActionChange;
    //cerr << "prevActionChange = " << prevActionChange << endl;
    //cerr << "logAcceptProb = " << logAcceptProb << endl;
    //cerr<<"My new action is "<<newAction<<" and my old action was
    //"<<oldAction<<endl;
    //    cout<<"The log of the accept prob is " << logAcceptProb;
    //    cout<<" "<<oldAction<<" "<<newAction<<" "<<endl;
    if (-logAcceptProb<log(sprng())){///reject conditin
      toAccept=false;
      //      break;

    }

    prevActionChange=currActionChange;
    levelCounter--;
  }
  if (toAccept ==true ){
    PathData.AcceptMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
    NumAccepted++;
    
    //cout<<"I'm accepting! I'm accepting!"<<endl;
  }
  else {
    PathData.RejectMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
    //  cout<<"I'm rejecting! I'm rejecting!"<<endl;
  }
  NumMoves++;
    

}
