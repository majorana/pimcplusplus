#include "BisectionClass.h"

/*! Samples the particles in "particles" and the time slices between
 startSlice and endSlice (not inclusive). Samples from a free gaussian
 of variance \f[$\sigma=\sqrt{2*\lambda*\tau} \$]. Return 
 the logarithm of the total sampling probability of the density (i.e \f[\$
\frac{-NDIM}{2}*\log{2*\pi*\sigma^2}-0.5*\frac{\delta \cdot
\delta}{\sigma * \sigma}  \$] )

*/
double BisectionClass::SamplePaths(int startSlice, int endSlice, Array<int,1> particles, int level)
{
  DistanceTableClass &distanceTable = *PathData.DistanceTable;
  dVec rpp;
  int skip = 1<<(level+1);
  double logNewSampleProb=0.0;
  PathClass &Path = PathData.Path;
  ActionClass &Action = PathData.Action;
  double levelTau = 0.5*PathData.Action.tau*skip;
  for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
    int  ptcl=particles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=startSlice;slice<endSlice;slice+=skip){
      dVec r = Path(slice,ptcl);
      dVec rp= Path(slice+skip,ptcl);
      dVec rpp=Path(slice+(skip>>1),ptcl);
      dVec rdiff=distanceTable.Velocity(slice, slice+skip, ptcl);
      //dVec rdiff = Path(slice+skip,ptcl)-Path(slice,ptcl);
      dVec rbar = r + 0.5*rdiff;
      dVec newDelta;
      Path.Random.LocalGaussianVec(sigma,newDelta);
      double GaussProd=1.0;
      for (int dim=0; dim<NDIM; dim++) {
	  while (newDelta[dim] > (0.5*Path.Box[dim]))
	    newDelta[dim] -= Path.Box[dim];
	  while (newDelta[dim] < (-(0.5*Path.Box[dim])))
	    newDelta[dim] += Path.Box[dim];
	  double GaussSum = 0.0;
	  int NumImage = 4;
	  for (int image=-NumImage; image <= NumImage; image++) {
	    double dist = newDelta[dim]+(double)image*Path.Box[dim];
	    GaussSum += exp(-0.5*dist*dist/sigma2);
	  }
	  GaussProd *= GaussSum;
      }
      logNewSampleProb += prefactorOfSampleProb + log(GaussProd);
      //distanceTable.PutInBox(newDelta);
      rpp=rbar+newDelta;
      //logNewSampleProb=logNewSampleProb+
      //	(prefactorOfSampleProb-0.5*dot(newDelta,newDelta)/(sigma2));
      ///Here we've stored the new position in the path
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
    }
  }
  return logNewSampleProb;
}


/// This calculates the sample probability for going from the state
/// that is currently in the newMode of MirroredArrayClass to the
/// state that is currently in oldMode of MirroredArrayClass 
double BisectionClass::LogSampleProb(int startSlice, int endSlice, 
				     Array<int,1> particles, int level)
{
  DistanceTableClass &distanceTable = *PathData.DistanceTable;
  PathClass &Path = PathData.Path;
  dVec rpp;
  int skip = 1<<(level+1);
  double logSampleProb=0.0;

  double levelTau = 0.5*PathData.Action.tau*skip;
  for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
    int ptcl = particles(ptclIndex);
    double lambda=Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=startSlice;slice<endSlice;slice+=skip){
      dVec r = Path(slice,ptcl);
      dVec rdiff = distanceTable.Velocity(slice, slice+skip, ptcl);
      //dVec rdiff = Path(slice+skip,ptcl)-Path(slice,ptcl);
      dVec rp= Path(slice+skip,ptcl);
      dVec rpp=Path(slice+(skip>>1),ptcl);
      ///We've ignored boundary conditions here (well we think this is fixed but we're not sure)
      dVec rbar=r + 0.5*rdiff;
      dVec Delta= rpp - rbar;
      distanceTable.PutInBox(Delta);
      
      double GaussProd=1.0;
      for (int dim=0; dim<NDIM; dim++) {
	double GaussSum = 0.0;
	int NumImage = 4;
	for (int image=-NumImage; image <= NumImage; image++) {
	  double dist = Delta[dim]+(double)image*Path.Box[dim];
	  GaussSum += exp(-0.5*dist*dist/sigma2);
	}
	GaussProd *= GaussSum;
      }
      logSampleProb += prefactorOfSampleProb + log(GaussProd);

      //logSampleProb=logSampleProb+
      //	(prefactorOfSampleProb-0.5*dot(Delta,Delta)/(sigma2));
    }
  }

  return logSampleProb;
}



// This does a multilevel construction of a path.
bool BisectionClass::Bisect(int startSlice,int numLevels, Array<int,1> activeParticles)
{
  bool toAccept=true;
  double oldLogSampleProb;
  double newLogSampleProb;
  Array<int,1> theParticles;
  int endSlice=(1<<numLevels)+startSlice;
  double prevActionChange=0;

  int levelCounter=numLevels-1;

  while (levelCounter>=0 && toAccept==true){
    SetMode(OLDMODE);
    toAccept=true;
    double oldAction = PathData.Action.calcTotalAction
      (startSlice,endSlice,activeParticles,levelCounter);
    oldLogSampleProb = LogSampleProb
      (startSlice,endSlice,activeParticles,levelCounter);
    SetMode(NEWMODE);
    newLogSampleProb = 
      SamplePaths(startSlice,endSlice,activeParticles,levelCounter);
    double testNewLogSampleProb=
      LogSampleProb(startSlice,endSlice,activeParticles,levelCounter);
    PathData.Update(startSlice,endSlice,activeParticles,
		    levelCounter);

    double newAction = PathData.Action.calcTotalAction
      (startSlice,endSlice, activeParticles,levelCounter);
    double currActionChange=newAction-oldAction;
    double logAcceptProb=
      -oldLogSampleProb+newLogSampleProb+currActionChange-prevActionChange;

    if (-logAcceptProb<log(PathData.Path.Random.Local())) ///reject conditin
      toAccept=false;

    prevActionChange=currActionChange;
    levelCounter--;
  }
  return (toAccept);
}
