#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "Common.h"
#include "SpeciesClass.h"


void PrintMoveClass::Read(IOSectionClass &IO)
{
  string typeCheck;
  assert(IO.ReadVar("type",typeCheck));
  assert(typeCheck=="PrintMove");
  assert(IO.ReadVar("name",Name));
  assert(IO.ReadVar("toprint",MyString));
}



void ShiftMoveClass::Read(IOSectionClass &theInput)
{
  string typeCheck;
  assert(theInput.ReadVar("type",typeCheck));
  assert(typeCheck=="ShiftMove");
  assert(theInput.ReadVar("name",Name));

}


void ShiftMoveClass::MakeMove()
{//Remember to mark Actions dirty!!!
  //int numTimeSlicesToShift=(int)floor(sprng()*PathData->NumTimeSlices);
  int numTimeSlicesToShift =(int)floor((PathData.NumTimeSlices()-1)*PathData.Path.Random.Common());
  
  //      PathData.MoveJoin(1);
  PathData.ShiftData(numTimeSlicesToShift);
  //      PathData.Join=1+numTimeSlicesToShift;

  
    
}

void BisectionMoveClass::Read(IOSectionClass &moveInput)
{
  string typeCheck;
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="Bisection");
  assert(moveInput.ReadVar("name",Name));
  Array<int,1> tempActiveSpecies(0);
  assert(moveInput.ReadVar("ActiveSpecies",tempActiveSpecies));
  SetActiveSpecies(tempActiveSpecies);
  int tempNumParticlesToMove;
  assert(moveInput.ReadVar("NumParticlesToMove",tempNumParticlesToMove));
  SetNumParticlesToMove(tempNumParticlesToMove);
  assert(moveInput.ReadVar("StartTimeSlice",StartTimeSlice));
  string tempNumLevels="";
  assert(moveInput.ReadVar("NumLevels",tempNumLevels));
  if (tempNumLevels=="Max"){
    NumLevels=PathData.Action.MaxLevels;
  }
  else {
    cerr<<"Don't know how to different number of levels from max yet"<<endl;
    assert(1==2);
  }

  //  moveInput->ReadVar("NumLevels",NumLevels);
  // int tempNumParticlesToMove;
  //moveInput->ReadVar("NumParticlesToMove",tempNumParticlesToMove);
  //SetNumParticlesToMove(tempNumParticlesToMove);
  //StartTimeSlice=0;

  ///HACK! HACK! HACK! HACK! Have to Find right way to input this.
  //Array<int,1> activeSpecies(1);
  //activeSpecies(0) = 0;
  //SetActiveSpecies(activeSpecies);


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
    double testNewLogSampleProb= PathData.Action.LogSampleProb
      (StartTimeSlice,EndTimeSlice,ActiveParticles,levelCounter);
    PathData.Update(StartTimeSlice,EndTimeSlice,ActiveParticles,
		    levelCounter);
		      
    //cerr << "newLogSampleProb = " << newLogSampleProb << endl;
    //cerr << "oldLogSampleProb = " << oldLogSampleProb << endl;
    double newAction = PathData.Action.calcTotalAction
      (StartTimeSlice,EndTimeSlice, ActiveParticles,levelCounter);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=
      -oldLogSampleProb+newLogSampleProb+currActionChange-prevActionChange;

    //    cerr << "prevActionChange = " << prevActionChange << endl;
    //    cerr << "logAcceptProb = " << logAcceptProb << " "<<oldLogSampleProb<<endl;
    //    cerr<<"My new action is "<<newAction<<" and my old action was ";
    //    cerr<<oldAction<<endl;
    //cout<<"The log of the accept prob is " << logAcceptProb;
    //    cout<<" "<<oldAction<<" "<<newAction<<" "<<endl;
    if (-logAcceptProb<log(PathData.Path.Random.Local())){///reject conditin
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
    //    cout<<"I'm rejecting! I'm rejecting!"<<endl;
  }
  NumMoves++;
    

}
