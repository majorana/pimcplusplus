#include "PathDataClass.h"
#include "BisectionMoveClass2.h"
#include "Common.h"
#include "SpeciesClass.h"


void BisectionMoveClass2::Read(IOSectionClass &moveInput)
{
  string typeCheck;
  assert(moveInput.ReadVar("type",typeCheck));
  assert(typeCheck=="Bisection2");
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

}

    

void BisectionMoveClass2::MakeMove()
{
  int EndTimeSlice=(1<<NumLevels)+StartTimeSlice;
  ChooseParticles(); 
  bool toAccept=Bisection.Bisect(StartTimeSlice,NumLevels,ActiveParticles);
  if (toAccept ==true ){
    PathData.AcceptMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
    NumAccepted++;
  }
  else {
    PathData.RejectMove(StartTimeSlice,EndTimeSlice,ActiveParticles);
  }
  NumMoves++;
}
