#include "HermeleFourierMove.h"
//HACK! HACK! HACK! THere is a 450 in this code which should be NumTimeSlices()-1
double HermeleFourierStageClass::Sample (int &slice1, int &slice2,
				   Array<int,1> &activeParticles)
{
#ifdef ORDER_N_FERMIONS


  double T=1.0/(PathData.Path.tau*(PathData.Path.NumTimeSlices()-1));
  PathData.Path.Phi2Omega();
  
  int randomPhi=PathData.Path.Random.LocalInt((PathData.Path.NumTimeSlices()-1)/2);
  double moveOmega=PathData.Path.Random.Local()-0.5;
  double moveImaginaryOmega=PathData.Path.Random.Local()-0.5;
  PathData.Path.outOmega[randomPhi][0]+=moveOmega;
  if (randomPhi!=0){
  PathData.Path.outOmega[450-randomPhi][0]+=moveOmega;
    PathData.Path.outOmega[randomPhi][1]+=moveImaginaryOmega;
    PathData.Path.outOmega[450-randomPhi][1]-=moveImaginaryOmega;
  }
  PathData.Path.Omega2Phi();
  for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++){
    PathData.Path(slice,0)[0]=PathData.Path.outPhi[slice][0];
    if (abs(PathData.Path.outPhi[slice][1])>1e-10){
      cerr<<PathData.Path.outPhi[slice][0]<<" "<<PathData.Path.outPhi[slice][1]<<endl;
      cerr<<randomPhi<<endl;
      cerr<<moveOmega<<" "<<moveImaginaryOmega<<endl;
    }

    assert(abs(PathData.Path.outPhi[slice][1])<1e-10);
  }
    PathData.Path(PathData.NumTimeSlices()-1,0)=PathData.Path(0,0);
  // And return sample probability ratio
#endif
  return 1.0;
}

void
HermeleFourierMoveClass::Read (IOSectionClass &in)
{
  //  cerr<<"Calling read of variational displace move"<<endl;
  // Construct action list
  //  assert(in.ReadVar("Sigma",HermeleFourierStage.Sigma));
  //  HermeleFourierStage.Actions.push_back(&PathData.Actions.ShortRangeOn);
  //  HermeleFourierStage.Actions.push_back(&PathData.Actions.VariationalPI);
  HermeleFourierStage.Actions.push_back(&PathData.Actions.DualHermele);


  
  // Now construct stage list
  Stages.push_back(&HermeleFourierStage);

  ActiveParticles.resize(1);

}



void
HermeleFourierMoveClass::MakeMove ()
{
  
  // First, choose particle to move
  int chosenPtcl =PathData.Path.Random.CommonInt(PathData.Path.NumParticles());
  ActiveParticles(0) = chosenPtcl;
  assert(chosenPtcl<PathData.Path.NumParticles());
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;
  //  cerr<<"Move begun"<<endl;
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;

  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }
  //  cerr<<"After levi flight"<<endl;
  //  PathData.Path.PrintRealSlices();

  TimesCalled++;
  if (toAccept)
    MultiStageClass::Accept();
  else 
    MultiStageClass::Reject();
  //  cerr<<"Move done"<<endl;
  //MoveClass::MakeMove();
}
