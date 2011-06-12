#include "WormSwap.h"

void
WormSwapMoveClass::Read (IOSectionClass &in)
{
  // Construct action list
  // WormSwapStage.Actions.push_back(&PathData.Actions.ShortRange);
  
  // Now construct stage list
  WormSwapStage.Read(in);
  WormSwapStage.Actions.push_back(&PathData.Actions.Kinetic);
  
  Stages.push_back(&WormSwapStage);

  ActiveParticles.resize(1);
}



///Warning: Moves Join
void 
WormSwapMoveClass::MakeRoomForGrowingTail(int slicesToGrow,int &tailSlice)
{
  int lastSlice=PathData.Path.NumTimeSlices()-1;
  PathData.MoveJoin(lastSlice);
  if (lastSlice-tailSlice>slicesToGrow)
    return;
  int shiftNeeded=slicesToGrow+tailSlice-lastSlice+1;
  PathData.Path.ShiftData(-shiftNeeded);
  PathData.Join=lastSlice-shiftNeeded;
  PathData.MoveJoin(lastSlice);
  tailSlice=tailSlice-shiftNeeded;

}

void 
WormSwapMoveClass::ChooseTimeSlices()
{
  int bisectSize = 1<<NumLevels;
  int tailSlice,tailPtcl;
  PathData.FindTail(tailSlice,tailPtcl);
  MakeRoomForGrowingTail(bisectSize,tailSlice);
  PathData.FindTail(tailSlice,tailPtcl);
  Slice1=tailSlice;
  Slice2=tailSlice+bisectSize;
}

int
WormSwapMoveClass::ChooseParticles()
{
  int swapPtcl=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
  //This works because you have empty particles after the tail and so
  //it's never the case that you can "pick" the head
  while (PathData.Path.ParticleExist(Slice2,swapPtcl)==0) 
    swapPtcl=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
  return swapPtcl;


}



void WormSwapMoveClass::MakeMove()
{
  int tailSlice,tailPtcl;
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);
  int permutePtcl=ChooseParticles();
  PathData.FindTail(tailSlice,tailPtcl);
  dVec tempVec=PathData.Path(Slice2,tailPtcl);
  PathData.Path(Slice2,tailPtcl)=PathData.Path(Slice2,permutePtcl);
  PathData.Path(Slice2,permutePtcl)=tempVec;
  ActiveParticles.resize(2);
  ActiveParticles(0)=

//   if (PathData.Path.OrderN){
//     for (int slice=Slice1;slice<=Slice2;slice++)
//       PathData.Path.Cell.BinParticles(slice);
//   }

  ((PermuteStageClass*)PermuteStage)->InitBlock(Slice1,Slice2);
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0)=-1;
    MultiStageClass::MakeMove();
  }

  if (LowestLevel != 0)
    MakeStraightPaths();

}
