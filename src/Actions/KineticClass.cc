#include "KineticClass.h"

///This has to be called after pathdata knows how many
///particles it has
void KineticClass::Read(IOSectionClass& in)
{
  cerr<<"I'm about to resize things now"<<endl;
  DoPtcl.resize(PathData.Path.NumParticles());
  cerr<<"I've finished resizing things"<<endl;
}

KineticClass::KineticClass(PathDataClass &pathData,
				 Array<PairActionFitClass* ,2> &pairMatrix) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix)
{
}

double KineticClass::Action (int slice1, int slice2,
				const Array<int,1> &changedParticles,
				int level)
{
  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      dVec vel;
      vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
      double GaussProd = 1.0;
      for (int dim=0; dim<NDIM; dim++) {
	int NumImage=1;
	double GaussSum=0.0;
	for (int image=-NumImage; image<=NumImage; image++) {
	  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	  GaussSum += exp(-dist*dist*FourLambdaTauInv);
	}
	GaussProd *= GaussSum;
      }
      TotalK -= log(GaussProd);    
      //TotalK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
  return (TotalK);
}



double KineticClass::d_dBeta (int slice1, int slice2,
				 int level)
{
  double spring=0.0;
  double levelTau=Path.tau;
  int slice2 = slice1 + (1<<level);
  for (int i=0; i<level; i++) 
    levelTau *= 2.0;
  spring  = 0.0;
  const int NumImage=1;
  for (int ptcl=0; ptcl<numPtcls; ptcl++)
    if (PathData.Path.ParticleSpecies(ptcl).lambda != 0.0)
      spring += 1.5/levelTau;
  
  for (int ptcl=0; ptcl<numPtcls; ptcl++) {
    // Do free-particle part
    int species1 = PathData.Path.ParticleSpeciesNum(ptcl);
    double lambda = PathData.Path.ParticleSpecies(ptcl).lambda;
    if (lambda != 0.0) {
      double FourLambdaTauInv = 
	1.0/(4.0*PathData.Path.Species(species1).lambda*levelTau);
      dVec vel;
      vel = PathData.Path.Velocity(slice1, slice2, ptcl);
      double Z = 1.0;
      dVec GaussSum=0.0;
      for (int dim=0; dim<NDIM; dim++) {
	for (int image=-NumImage; image<=NumImage; image++) {
	  double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	  GaussSum[dim] += exp(-dist*dist*FourLambdaTauInv);
	}
	Z *= GaussSum[dim];
      }
      dVec numSum=0.0;
      for (int dim=0;dim<NDIM;dim++){
	for (int image=-NumImage;image<=NumImage;image++){
	  double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	  numSum[dim] += 
	    (-dist*dist*FourLambdaTauInv/levelTau)*exp(-dist*dist*FourLambdaTauInv);
	}
      }
      double scalarnumSum=0.0;
      for (int dim=0;dim<NDIM;dim++){
	dVec numProd=1.0;
	for (int dim2=0;dim2<NDIM;dim2++){
	  if (dim2!=dim){
	    numProd[dim] *= GaussSum[dim2];
	  }
	  else {
	    numProd[dim] *=  numSum[dim2];
	  }
	  
	}
	scalarnumSum += numProd[dim];
      }
      spring += scalarnumSum/Z; 
    }
  }
    
  return spring;
}
