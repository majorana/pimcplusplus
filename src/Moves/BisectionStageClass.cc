#include "BisectionStageClass.h"

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

