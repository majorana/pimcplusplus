#include "KineticClass.h"
#include "../PathDataClass.h"

///This has to be called after pathdata knows how many
///particles it has
void KineticClass::Read(IOSectionClass& in)
{
}

KineticClass::KineticClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double KineticClass::Action (int slice1, int slice2,
			     const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  Array<double,1> KineticVal(PathData.Path.NumTimeSlices());
  KineticVal=0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double lambda = Path.Species(species).lambda;
    if (lambda != 0){
      double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
      for (int slice=slice1; slice < slice2;slice+=skip) {
        dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	//vel= PathData.Path(slice+skip,ptcl)-PathData.Path(slice,ptcl);
	
        double GaussProd = 1.0;
        for (int dim=0; dim<NDIM; dim++) {
	  double GaussSum=0.0;
	  for (int image=-NumImages; image<=NumImages; image++) {
	    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	    GaussSum += exp(-dist*dist*FourLambdaTauInv);
	  }
	  GaussProd *= GaussSum;
        }
	TotalK -= log(GaussProd);    
	KineticVal(slice)-=log(GaussProd);
	//TotalK += dot(vel,vel)*FourLambdaTauInv; 
	//	KineticVal(slice)+=dot(vel,vel)*FourLambdaTauInv; 
      }
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
  //  cerr << "I'm returning kinetic action " << TotalK << endl;
  //  for (int counter=0;counter<KineticVal.size();counter++){
  //    cerr<<"My kinetic link "<<counter<<" is "<<KineticVal(counter)<<endl;
  //  }
  return (TotalK);
}



double KineticClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  Array<double,1> KineticVal(PathData.Path.NumTimeSlices());
  KineticVal=0.0;

  double spring=0.0;
  // ldexp(double x, int n) = x*2^n
  double levelTau=ldexp(Path.tau, level);
//   for (int i=0; i<level; i++) 
//     levelTau *= 2.0;
  spring  = 0.0;  
  int skip = 1<<level;
  for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++) {
    // Do free-particle part
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    SpeciesClass &species = Path.Species(speciesNum);
    double lambda = species.lambda;
    if (lambda != 0.0) {
      double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice<slice2; slice+=skip) {
	spring += (0.5*NDIM)/levelTau;
	KineticVal(slice)+=(0.5*NDIM)/levelTau;
	dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	double Z = 1.0;
	dVec GaussSum=0.0;
	dVec numSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImages; image<=NumImages; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    double d2overFLT = dist*dist*FourLambdaTauInv;
	    double expPart = exp(-d2overFLT);
	    GaussSum[dim] += expPart;
	    numSum[dim] += -d2overFLT/levelTau* expPart;
	  }
	  Z *= GaussSum[dim];
	}
	double scalarnumSum=0.0;
	for (int dim=0;dim<NDIM;dim++) {
	  dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
	} //cerr << "Z = " << Z << " scalarnumSum = " << scalarnumSum << endl;
	spring += scalarnumSum/Z; 
	KineticVal(slice)+=scalarnumSum/Z;
      }
    }
  }
  //  cerr << "spring = " << spring << endl;
  //cerr << "I'm returning kinetic energy " << TotalK << endl;
//   cerr<<"Returning kinetic energy"<<endl;
//   for (int counter=0;counter<KineticVal.size();counter++){
//     cerr<<"My kinetic link "<<counter<<" is "<<KineticVal(counter)<<endl;
//   }

  return spring;
}
