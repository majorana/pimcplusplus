#include "WormGrow.h"

void
WormGrowMoveClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar("MaxToGrow",MaxGrowth));
  // Construct action list
  //  WormGrowStage.Actions.push_back(&PathData.Actions.ShortRange);
  
  // Now construct stage list
  WormGrowStage.Read(in);
  WormGrowStage.Actions.push_back(&PathData.Actions.Kinetic);
  WormGrowStage.Actions.push_back(&PathData.Actions.Mu);
  
  Stages.push_back(&WormGrowStage);
  
  ActiveParticles.resize(1);
}

void 
WormGrowMoveClass::MakeMove()
{
  Slice1=0;
  Slice2=0;
  ActiveParticles.resize(1);
  ActiveParticles(0)=0;
  cerr<<"Move begun"<<endl;
  if (!PathData.Path.NowOpen){
    MultiStageClass::Reject();
    return;
  }
  int headSlice,headPtcl,tailSlice,tailPtcl,numEmpty,wormSize;
  PathData.WormInfo(headSlice, headPtcl,
		    tailSlice, tailPtcl,
		    numEmpty, wormSize);
  
  cerr<<"Worm size is "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;
  if (wormSize==0){
    MultiStageClass::Reject();
    return;
  }
  bool moveHead = PathData.Path.Random.Local() >0.5;
  bool grow = PathData.Path.Random.Local() < 0.5;
  int changeAmount=PathData.Path.Random.LocalInt(MaxGrowth)+1;
  //  cerr<<"Change amount is "<<changeAmount<<endl;
  if (grow && numEmpty-changeAmount<PathData.Path.NumTimeSlices()){
    MultiStageClass::Reject();
    return;
  }
  
  if (!grow && wormSize-changeAmount<1){
    MultiStageClass::Reject();
    return;
  }
  WormGrowStage.MoveHead=moveHead;
  WormGrowStage.Grow=grow;
  WormGrowStage.ChangeAmount=changeAmount;

  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;

  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }


  PathData.WormInfo(headSlice, headPtcl,
		    tailSlice, tailPtcl,
		    numEmpty, wormSize);
  
  cerr<<"Worm size as "<<wormSize<<" "<<numEmpty<<" "
      <<headSlice<<" "<<headPtcl<<" "
      <<tailSlice<<" "<<tailPtcl<<" "<<endl;

//   cerr<<"Worm information: "<<headSlice<<" "<<headPtcl<<" "
//       <<tailSlice<<" "<<tailPtcl<<" "
//       <<numEmpty<<" "<<wormSize<<endl;
  assert(numEmpty>=PathData.Path.NumTimeSlices());
  assert(wormSize>0);
  

  TimesCalled++;
  if (toAccept)
    Accept();
  else 
    Reject();

  
  cerr<<"Move done"<<endl;

  

}
