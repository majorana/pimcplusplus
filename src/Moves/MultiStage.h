#ifndef MULTI_STAGE_H
#define MULTI_STAGE_H


#include "list.h"
#include "MoveClass.h"

class ActionBase
{
public:
  double Action();

};

class StageClass
{
protected:
  PathDataClass &PathData;
public:
  list<ActionBase*> Actions;
  ///The highest stage will set the slices and activeParticles
  ///This returns transition probability T(new->old)/T(old->new)
  virtual double Sample (int &slice1,int &slice2,
			 Array<int,1> activeParticles); 
  double StageAction(int startSlice,int endSlice,
		     const Array<int,1> &changedParticles);
  StageClass(PathDataClass &pathData) :PathData(pathData)
  {
    //do nothing for now
  }
};


double StageClass::StageAction(int startSlice,int endSlice,
		   const Array<int,1> &changedParticles)
{
  double TotalAction=0.0;
  list<ActionBase*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()){
    TotalAction += ((*actionIter)->Action());
    actionIter++;
  }
  return TotalAction;

}


class BisectionStageClass : public StageClass
{
private:
  int Level;
  int NumImage;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> activeParticles);
  BisectionStageClass(PathDataClass &pathData,int level) : StageClass(pathData),
							   Level(level)
  { 
    //do nothing for now
  }
  



};

double BisectionStageClass::Sample(int &slice1,int &slice2,
				   Array<int,1> activeParticles)
{
  PathClass &Path = PathData.Path;
  ActionClass &Action = PathData.Action;
  int skip = 1<<(Level+1);
  double levelTau = 0.5*PathData.Action.tau*skip;

  double logSampleProb=0.0;
  double logOldSampleProb=0.0;

  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=slice1;slice<slice2;slice+=skip){
      SetMode(OLDMODE);
      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+ 0.5*rdiffOld;
      dVec rppOld=Path(slice+(skip>>1),ptcl);
      dVec DeltaOld=rppOld-rbarOld;
      Path.PutInBox(DeltaOld);
      
      SetMode(NEWMODE);
      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbar=r+ 0.5*rdiff;
      dVec  Delta;
      Path.Random.LocalGaussianVec(sigma,Delta);
      PathData.Path.PutInBox(Delta);

      double GaussProd=1.0;
      double GuassProdOld=1.0;
      for (int dim=0; dim<NDIM; dim++) {
	double GaussSum = 0.0;
	double GaussSumOld =0.0;
	for (int image=-NumImage; image <= NumImage; image++) {
	  double dist = Delta[dim]+(double)image*Path.GetBox()[dim];
	  double distOld=DeltaOld[dim]+(double)image*Path.GetBox()[dim];
	  GaussSum += exp(-0.5*dist*dist/sigma2);
	  GaussSumOld += exp(-0.5*distOld*distOld/sigma2);
	}
	GaussProd *= GaussSum;
	GuassProdOld *= GaussSumOld;
      }
      logOldSampleProb += prefactorOfSampleProb + log(GaussProdOld);
      logSampleProb += prefactorOfSampleProb  + log(GaussProd);
    
    rpp=rbar+newDelta;
    ///Here we've stored the new position in the path
    Path.SetPos(slice+(skip>>1),ptcl,rpp);
    }
  }
  return logSampleProb-logOldSampleProb;
}


class MultiStageClass : public ParticleMoveClass
{
protected:
  list<StageClass*> Stages;
  int NumSteps;
public:
  void Read(IOSectionClass &io);
  void MakeMove();
  MultiStageClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    ParticleMoveClass(pathData,outSection) 
  {
    //do nothing for now
  }
};

void MultiStageClass::MakeMove()
{

  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end()){
    SetMode(OLDMODE);
    double oldAction=stageIter->StageAction(slice1,slice2,activeParticles);
    SetMode(NEWMODE);
    double sampleRatio=stageIter->Sample(slice1,slice2,activeParticles);    
    double newAction = stageIter->StageAction(slice1,slice2,activeParticles)
    double currActionChange=newAction-oldAction;
    double logAcceptProb=sampleRatio+currActionChange-prevActionChange;
    if (-logAcceptProb<log(PathData.Path.Random.Local())) ///reject conditin
      toAccept=false;
    if (toAccept)
      NumAccepted(levelCounter)++;
    else
      NumRejected(levelCounter)++;
    prevActionChange=currActionChange;
    stageIter++;
  }
}

class 

#endif
