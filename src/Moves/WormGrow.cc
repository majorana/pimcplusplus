#include "WormGrow.h"

void
WormGrowMoveClass::Read (IOSectionClass &in)
{
  // Construct action list
  // WormGrowStage.Actions.push_back(&PathData.Actions.ShortRange);
  
  // Now construct stage list
  Stages.push_back(&WormGrowStage);

  ActiveParticles.resize(1);
}

double
WormGrowStageClass::Sample (int &slice1, int &slice2,
			    Array <int,1> &activeParticles)
{
  cerr<<"Sampling now"<<endl;
  // Move the Join out of the way.
  PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);
  cerr<<"IO've moved the join"<<endl;
  //only do helium
  int ptcl=0;
  cerr<<"Resizing"<<endl;
  activeParticles.resize(1);

  //  bool doHead=(PathData.Path.Random.Local()>0.5);
  bool doHead=true;
  bool doGrow=(PathData.Path.Random.Local()>0.5);
  int numSlices=PathData.Path.Random.LocalInt(MaxSlicesToGrow)+1;
  double logSampleProb=0.0;
  cerr<<"here"<<endl;
  if (doHead && doGrow){
    cerr<<"Do head and do grow"<<endl;
    while (PathData.Path.OpenLink+numSlices>=PathData.Path.NumTimeSlices()){
      PathData.MoveJoin(1);
      PathData.Path.ShiftData(-1);
      PathData.Join=0;
      PathData.Path.HeadSlice=
	(PathData.Path.HeadSlice-1+PathData.Path.NumTimeSlices()) % 
	PathData.Path.NumTimeSlices();
    }
    cerr<<"Post do head stuff"<<endl;
    PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);
    cerr<<"Moved join"<<endl;
    activeParticles(0)=PathData.Path.HeadParticle;    
    double tau=PathData.Path.tau;
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(2.0*lambda*tau);
    double sigma=sqrt(sigma2);
    dVec delta;
    cerr<<"hmm"<<endl;
    for (int counter=0;counter<numSlices;counter++){
      cerr<<"A"<<endl;
      cerr<<PathData.Path.HeadSlice<<" "<<PathData.Path.HeadParticle<<endl;
      
      PathData.Path.Random.LocalGaussianVec(sigma,delta);
      //      PathData.Path.PushHead(PathData.Path.CurrHead+delta);
      PathData.Path(PathData.Path.HeadSlice+1,PathData.Path.HeadParticle)=
	PathData.Path(PathData.Path.HeadSlice,PathData.Path.HeadParticle)+delta;
      PathData.Path.HeadSlice=(PathData.Path.HeadSlice+1) %
	PathData.Path.NumTimeSlices();
      logSampleProb+=(-0.5*delta*delta/sigma2);
    }
    cerr<<"end a"<<endl;
    slice1 = 0;
    slice2 = PathData.Path.NumTimeSlices()-1;
    cerr<<"Getting out once"<<endl;
    return exp(logSampleProb);
  }
  else {
    cerr<<"Chose this path"<<endl;
    while (PathData.Path.OpenLink-numSlices>=0){
      PathData.MoveJoin(0);
      PathData.Path.ShiftData(1);
      PathData.Join=1;
      PathData.Path.HeadSlice=
	(PathData.Path.HeadSlice+1+PathData.Path.NumTimeSlices()) % 
	PathData.Path.NumTimeSlices();

    }
    cerr<<"Going to move join"<<endl;
    PathData.MoveJoin (0);
    
    activeParticles(0)=PathData.Path.HeadParticle;    
    for (int counter=0;counter<numSlices;counter++){
      PathData.Path.HeadSlice=PathData.Path.HeadSlice-1;
    }
    cerr<<"Helooo"<<endl;
    slice1 = 0;
    slice2 = PathData.Path.NumTimeSlices()-1;
    cerr<<"Getting out twice"<<endl;
    return 1.0;


  }
  

}
