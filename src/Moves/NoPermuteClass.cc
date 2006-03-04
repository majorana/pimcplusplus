#include "NoPermuteStage.h"

double NoPermuteStageClass::Sample (int &slice1, int &slice2,
				    Array<int,1> &activeParticles)
{
  return 1.0;
}

bool NoPermuteStageClass::Attempt (int &slice1, int &slice2,
				   Array<int,1> &activeParticles, double &prevActionChange)
{
  // If I'm being called the first time in the move, choose a particle
  if (activeParticles(0) == -1) {
    activeParticles.resize (1);
    //HACK! HACK! HACK! HACK!
    if (PathData.Path.OpenPaths && slice1<PathData.Path.OpenLink &&
	PathData.Path.OpenLink<slice2 && 1==2){
      do {
	activeParticles(0) = ChooseParticle();
      } while (activeParticles(0)==PathData.Path.OpenPtcl);
    }
    else activeParticles(0)=ChooseParticle();
 
  }
  // In any case, always accept
  return true;
}

void NoPermuteStageClass::InitBlock(int &slice1,int &slice2)
{
  //  cerr<<"In init block"<<endl;
  if (PathData.Path.OrderN){
    for (int slice=slice1;slice<=slice2;slice++)
      PathData.Path.Cell.BinParticles(slice);
  }
  //  cerr<<"out of init block"<<endl;

}
int NoPermuteStageClass::ChooseParticle()
{
  int myPtcl=(PathData.Path.Random.LocalInt 
    (PathData.Path.Species(SpeciesNum).NumParticles))+
    PathData.Path.Species(SpeciesNum).FirstPtcl;
  //  cerr<<"I've chosen the particle "<<myPtcl;
  return myPtcl;
}

