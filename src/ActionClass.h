#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "Common/Splines/CubicSpline.h"
#include "Common/PairAction/PAFit.h"
#include "SpeciesClass.h"
#include "PathClass.h"
//#include "DistanceTablePBCClass.h"
//#include "DistanceTableFreeClass.h"
#include "DistanceTableClass.h"



/// This is the class that controls all of the actions and is in
/// charge of calculating them. When this is initialized a pointer needs
/// to be sent that has the memoizedData and SpeciesClass 
class ActionClass
{
private:
public:
  DistanceTableClass *DistanceTable;
  /// This holds all of the Pair Action Classes
  Array<PairActionFitClass*,1> PairActionVector;
  /// Holds indices to which PairActionClass in the PairAcctionVector
  /// you use for a given pair of particles indexed by
  /// (species1,species2) 
  Array<int,2> PairMatrix;
  PathClass &Path;
  /// Temperature
  double tau;
  /// The maximum number of levels we can have in a bisection move.
  int MaxLevels;


  void Read(IOSectionClass &IOSection);
  ActionClass(PathClass  &p_path) : Path(p_path) 
  {
  }


  /// Calculates the total action.
  double calcTotalAction(int startSlice, int endSlice, 
			 Array<int,1> changedParticles,int level);

  ///This picks a new location in space for the particles in the
  ///particles Array at all of the time slices between startSlice and
  ///endSlice (at the appropriate skip for the level)

  inline double SampleParticles(int startSlice, int endSlice, 
				Array<int,1> particles, int level);

  /// This calculates the sample probability for going from the state
  /// that is currently in the newMode of MirroredArrayClass to the
  /// state that is currently in oldMode of MirroredArrayClass 
  inline double LogSampleProb(int startSlice, int endSlice, 
			      Array<int,1> particles,int level);

  /// Function to calculate the total action.
  void calcTotalAction();

};


/*! Samples the particles in "particles" and the time slices between
 startSlice and endSlice (not inclusive). Samples from a free gaussian
 of variance \f[$\sigma=\sqrt{2*\lambda*\tau} \$]. Return 
 the logarithm of the total sampling probability of the density (i.e \f[\$
\frac{-NDIM}{2}*\log{2*\pi*\sigma^2}-0.5*\frac{\delta \cdot
\delta}{\sigma * \sigma}  \$] )

*/
inline double ActionClass::SampleParticles(int startSlice, int endSlice, Array<int,1> particles, int level)
{
  dVec rpp;
  int skip = 1<<(level+1);
  double logNewSampleProb=0.0;

  double levelTau = 0.5*tau*skip;
  for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
    int  ptcl=particles(ptclIndex);
    double lambda=Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=startSlice;slice<endSlice;slice+=skip){
      dVec r = Path(slice,ptcl);
      dVec rp= Path(slice+skip,ptcl);
      dVec rpp=Path(slice+(skip>>1),ptcl);
      dVec rdiff=DistanceTable->Velocity(slice, slice+skip, ptcl);
      //dVec rdiff = Path(slice+skip,ptcl)-Path(slice,ptcl);
      dVec rbar = r + 0.5*rdiff;
      dVec newDelta=GaussianRandomVec(sigma);
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
      //DistanceTable->PutInBox(newDelta);
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
inline double ActionClass::LogSampleProb(int startSlice, int endSlice, 
					 Array<int,1> particles, 
					 int level)
{
  dVec rpp;
  int skip = 1<<(level+1);
  double logSampleProb=0.0;

  double levelTau = 0.5*tau*skip;
  for (int ptclIndex=0;ptclIndex<particles.size();ptclIndex++){
    int ptcl = particles(ptclIndex);
    double lambda=Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int slice=startSlice;slice<endSlice;slice+=skip){
      dVec r = Path(slice,ptcl);
      dVec rdiff = DistanceTable->Velocity(slice, slice+skip, ptcl);
      //dVec rdiff = Path(slice+skip,ptcl)-Path(slice,ptcl);
      dVec rp= Path(slice+skip,ptcl);
      dVec rpp=Path(slice+(skip>>1),ptcl);
      ///We've ignored boundary conditions here (well we think this is fixed but we're not sure)
      dVec rbar=r + 0.5*rdiff;
      dVec Delta= rpp - rbar;
      DistanceTable->PutInBox(Delta);
      
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



#endif
