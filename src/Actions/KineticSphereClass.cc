#include "KineticSphereClass.h"
#include "../PathDataClass.h"
#include "../Common/SpecialFunctions/SpecialFunctions.h"
///This has to be called after pathdata knows how many
///particles it has
void KineticSphereClass::Read(IOSectionClass& in)
{
  assert(in.ReadVar("SphereRadius",SphereRadius));
}

KineticSphereClass::KineticSphereClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double KineticSphereClass::Action (int slice1, int slice2,
			     const Array<int,1> &changedParticles, int level)
{
  //  cerr<<"I am in the action"<<endl;
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
	vel= PathData.Path(slice+skip,ptcl)-PathData.Path(slice,ptcl);
	dVec r1=PathData.Path(slice,ptcl);
	dVec r2=PathData.Path(slice+skip,ptcl);
	double a=SphereRadius;
	double angle;
	if (dot(r1,r2)/(a*a)<=1)
	  angle= acos(dot(r1,r2)/(a*a));
	else 
	  angle=0;
	double vel2=angle*angle*a*a;
	double toAddAction=vel2*FourLambdaTauInv;
	////	double toAddAction2=log(K(slice,slice+skip,ptcl,level,lambda));
	//	if (level==0)
	//	  cerr<<toAddAction<<" "<<-toAddAction2+7.37<<" "<<vel*vel*FourLambdaTauInv
	//	      <<" "<<slice<<" "<<slice+skip<<" "<<level<<" "
	//	      <<" "<<toAddAction+toAddAction2<<" "
	//	      <<" "<<-toAddAction2/toAddAction<<endl;
	
	 
	///	TotalK -=toAddAction2; 
	TotalK += vel2*FourLambdaTauInv;


	
	
//         double GaussProd = 1.0;
//         for (int dim=0; dim<NDIM; dim++) {
// 	  double GaussSum=0.0;
// 	  for (int image=-NumImages; image<=NumImages; image++) {
// 	    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
// 	    GaussSum += exp(-dist*dist*FourLambdaTauInv);
// 	  }
// 	  GaussProd *= GaussSum;
// 	}
//	TotalK -= log(GaussProd);    
//	KineticVal(slice)-=log(GaussProd);
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
	spring += (0.5*(NDIM-1))/levelTau;
	KineticVal(slice)+=(0.5*(NDIM-1))/levelTau;
	///END HACK!
	dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	//	double vel2=dot(vel,vel);
	dVec r1=PathData.Path(slice,ptcl);
	dVec r2=PathData.Path(slice+skip,ptcl);
	double a=SphereRadius;
	double angle= acos(dot(r1,r2)/(a*a));
	double vel2=angle*angle*a*a;
	//	cerr<<"My angle is "<<angle<<" ";
	//	cerr<<dot(r1,r2)<<" ";
	//	cerr<<dot(r1,r2)/(a*a)<<" ";
	//	cerr<<acos(dot(r1,r2)/(a*a))<<" ";
	//	cerr<<(1==1)<<" "<<(dot(r1,r2)/(a*a)<=1)<<" ";
	//	cerr<<r1<<" "<<r2<<endl;
	///////	//	double oldEnergy=K2(slice,slice+skip,ptcl,level,lambda)*lambda/(SphereRadius*SphereRadius);
	//////	//	double newEnergy=(0.5*(NDIM-1))/levelTau-vel2*FourLambdaTauInv/levelTau;
	/////	//      double newestEnergy=0.5*(NDIM-1)/levelTau-2*vel2*FourLambdaTauInv/levelTau;
	//	cerr<<newEnergy<<" "<<oldEnergy<<" "<<newEnergy-oldEnergy<<" "<<newEnergy/oldEnergy<<" "<<newestEnergy<<" "<<newestEnergy/oldEnergy<<endl;
	spring += -vel2*FourLambdaTauInv/levelTau;
// 	double Z = 1.0;
// 	dVec GaussSum=0.0;
// 	dVec numSum=0.0;
// 	for (int dim=0; dim<NDIM; dim++) {
// 	  for (int image=-NumImages; image<=NumImages; image++) {
// 	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
// 	    double d2overFLT = dist*dist*FourLambdaTauInv;
// 	    double expPart = exp(-d2overFLT);
// 	    GaussSum[dim] += expPart;
// 	    numSum[dim] += -d2overFLT/levelTau* expPart;
// 	  }
// 	  Z *= GaussSum[dim];
// 	}
// 	double scalarnumSum=0.0;
// 	for (int dim=0;dim<NDIM;dim++) {
// 	  dVec numProd=1.0;
// 	  for (int dim2=0;dim2<NDIM;dim2++) {
// 	    if (dim2!=dim)
// 	      numProd[dim] *= GaussSum[dim2];
// 	    else 
// 	      numProd[dim] *=  numSum[dim2];
// 	  }
// 	  scalarnumSum += numProd[dim];
// 	} //cerr << "Z = " << Z << " scalarnumSum = " << scalarnumSum << endl;
// 	spring += scalarnumSum/Z; 
//	KineticVal(slice)+=scalarnumSum/Z;
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

double KineticSphereClass::d_dBeta_old (int slice1, int slice2,
				    int level)
{
  Array<double,1> KineticVal(PathData.Path.NumTimeSlices());
  KineticVal=0.0;
  //  cerr<<"IN here"<<endl;
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
 	spring += K2(slice,slice+skip,ptcl,level,lambda)*
	  lambda/(SphereRadius*SphereRadius);
	//	spring += scalarnumSum/Z; 
      }
    }
  }

  double oldSpring = d_dBeta_old(slice1,slice2,level);
  cerr<<spring<<" "<<oldSpring<<" "<<spring-oldSpring<<" "<<spring/oldSpring<<endl;
  return spring;
}


double KineticSphereClass::K(int slice,int nextSlice,int ptcl, int level,double lambda){

  double numerator=0,denominator=0;
  double num,dnum;
  //  cerr<<"In kineticsphereclass"<<endl;
  dVec mySlice=PathData.Path(slice,ptcl),
    myNextSlice=PathData.Path(nextSlice,ptcl);
  double cosGamma;
  double levelTau=ldexp(Path.tau, level);
  int maxl=500;
  double ans;
//   //  do {
  cosGamma=dot(myNextSlice,mySlice)/sqrt(dot(myNextSlice,myNextSlice)*dot(mySlice,mySlice));
  for (int l=0;l<=maxl;l++)
    {
      num=l*(l+1);
      dnum=(2*l+1)*Legendre(l,cosGamma)*
	exp(-levelTau*num*lambda/(SphereRadius*SphereRadius));
      num*=dnum;
	
      numerator+=num;
      denominator+=dnum;
    }
  ans=numerator/denominator;
//     //  }
// //  while(ans<0);
  
// //    return numerator;
//     return ans;
  return denominator;
//   //remember to multiply by hbar2/(2mradius) when u exit the function
  
}

double KineticSphereClass::K2(int slice,int nextSlice,int ptcl, int level,double lambda){

  double numerator=0,denominator=0;
  double num,dnum;
  //  cerr<<"In kineticsphereclass"<<endl;
  dVec mySlice=PathData.Path(slice,ptcl),
    myNextSlice=PathData.Path(nextSlice,ptcl);
  double cosGamma;
  double levelTau=ldexp(Path.tau, level);
  int maxl=500;
  double ans;
//   //  do {
  cosGamma=dot(myNextSlice,mySlice)/sqrt(dot(myNextSlice,myNextSlice)*dot(mySlice,mySlice));
  for (int l=0;l<=maxl;l++)
    {
      num=l*(l+1);
      dnum=(2*l+1)*Legendre(l,cosGamma)*
	exp(-levelTau*num*lambda/(SphereRadius*SphereRadius));
      num*=dnum;
	
      numerator+=num;
      denominator+=dnum;
    }
  ans=numerator/denominator;
//     //  }
// //  while(ans<0);
  
// //    return numerator;
//     return ans;
  return ans;
//   //remember to multiply by hbar2/(2mradius) when u exit the function
  
}
