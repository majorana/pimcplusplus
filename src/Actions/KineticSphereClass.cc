#include "KineticSphereClass.h"
#include "../PathDataClass.h"
#include "../Common/SpecialFunctions/SpecialFunctions.h"
///This has to be called after pathdata knows how many
///particles it has
void KineticSphereClass::Read(IOSectionClass& in)
{
}

KineticSphereClass::KineticSphereClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double KineticSphereClass::Action (int slice1, int slice2,
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



double KineticSphereClass::d_dBeta (int slice1, int slice2,
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
	////HACK FOR HELIUM SPHERE!
	//	spring += (0.5*(NDIM-1))/levelTau;
	spring += K(slice1,slice1+1,ptcl,level);
	  //	spring += scalarnumSum/Z; 
      }
    }
  }

  return spring;
}


double KineticSphereClass::K(int slice,int nextSlice,int ptcl, int level){

  double numerator=0,denominator=0;
  double num,dnum;
  dVec mySlice=PathData.Path(slice,ptcl),
    myNextSlice=PathData.Path(nextSlice,ptcl);
  double cosGamma;
  double levelTau=ldexp(Path.tau, level);
  int maxl=50;
  double ans;
  do {
    for (int l=0;l<=maxl;l++)
      {
	cosGamma=dot(myNextSlice,mySlice)/(dot(myNextSlice,myNextSlice)*dot(mySlice,mySlice));
	num=l*(l+1);
	dnum=(2*l+1)*Legendre(l,cosGamma)*exp(-levelTau*num);
	num*=dnum;
	
	numerator+=num;
	denominator+=dnum;
      }
    ans=numerator/denominator;
  }
  while(ans<0);
    
    return ans;
  
  //remember to multiply by hbar2/(2mradius) when u exit the function
  
}
