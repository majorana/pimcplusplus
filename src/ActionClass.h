#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "Common/Splines/CubicSpline.h"
#include "SpeciesClass.h"
#include "MemoizedDataClass.h"
#include "PathClass.h"
//#include "DistanceTablePBCClass.h"
//#include "DistanceTableFreeClass.h"
#include "DistanceTableClass.h"

// #include "SpeciesArrayClass.h"
// #include "ObservableClass.h"

/// Doesn't do anything right now but saves the pair action class
/// eventually 
class SavedPairActionClass
{
  ///This is a test comment


};

///This is the pair action class. It uses the following formula in
///order to calculate the pair action
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
class PairActionClass
{

 private:
  
  ///Holds the Ukj coefficients for a given q
  Array<double,1> TempukjArray;
  DistanceTableClass *DistanceTable;
  /// Skips to the next string in the file whose substring matches skipToString
  string SkipTo(ifstream &infile, string skipToString);
  /// Reads a Fortran 3 tensor
  void ReadFORTRAN3Tensor(ifstream &infile, Array<double,3> &tempUkj);
 public:
  string type1,type2;
  void Read(InputSectionClass &inSection);
  /// This stores the coefficients of the expansion specified above.
  /// The array index is over the levels.  You call it with the q
  /// value and a temporary array to get all of the values in that
  /// column. 
  Array<MultiCubicSpline,1> ukj; ///<(level )
  ///Same as ukj but stores the beta derivatives.
  Array<MultiCubicSpline,1> dukj; ///<(level )
  /// Calculate the U(s,q,z) value when given s,q,z and the level 
  inline double calcUsqz(double s,double q,double z,int level);
  /// This is the order of the fit to use. 
  int n;
  /// This is the temperature 
  double tau;
  /// Function to read David's squarer file input.
  void ReadDavidSquarerFile(string DMFile);

};

/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
inline double PairActionClass::calcUsqz(double s,double q,double z,int level)
{
  double sum=0.0;
  double r=q+0.5*z;
  double rprime=q-0.5*z;

  if (q > ukj(level).grid->End)
    return (0.0);


  if (r > ukj(level).grid->End)
    return (0.0);

  if (rprime > ukj(level).grid->End)
    return (0.0);



  sum=sum+0.5*((ukj(level))(0,r)+(ukj(level))(0,rprime)); //This is the endpoint action 


  if (s > 0.0)
    {
      double zsquared=z*z;
      double ssquared=s*s;
      double ssquaredinverse=1.0/ssquared;
      double Sto2k=ssquared;
      (ukj(level))(q,TempukjArray); 
      for (int k=1;k<=n;k++){  
	
	double Zto2j=1;
	double currS=Sto2k;
	
	for (int j=0;j<=k;j++){
	  
	  double cof=TempukjArray(k*(k+1)/2+j); //indexing into the 2darray
	  sum=sum+cof*Zto2j*currS;
	  
	  
	  Zto2j*=zsquared;
	  currS=currS*ssquaredinverse;				
	}				
	Sto2k=Sto2k*ssquared;
      }
    }
  
  return sum; 
}

/// This is the class that controls all of the actions and is in
/// charge of calculating them. When this is initialized a pointer needs
/// to be sent that has the memoizedData and SpeciesClass 
class ActionClass
{
private:
public:
  void Read(InputSectionClass &inSection);
  DistanceTableClass *DistanceTable;
  /// This holds all of the Pair Action Classes
  Array<PairActionClass,1> PairActionVector;
  /// Holds indices to which PairActionClass in the PairAcctionVector
  /// you use for a given pair of particles indexed by
  /// (species1,species2) 
  Array<int,2> PairMatrix;
  /// This indexes into the non-existent Saved Pair Actions
  Array<SavedPairActionClass,2> SavedPairActionArray;
  ActionClass(PathClass  &p_path) : Path(p_path) 
  {
  }
  /// This holds a reference to the Array of Species
  PathClass &Path;
  /// Temperature
  double tau;
  /// Calculates the total action.
  double calcTotalAction(int startSlice, int endSlice, 
			 Array<int,1> changedParticles,int level);
  /// This is a reference to the memoized data class
  //  MemoizedDataClass &myMemoizedData;
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
    //    int species=particles(ptcl)(0);
    //    int ptclNum=particles(ptcl)(1);
    double lambda=Path.ParticleSpecies(ptcl).lambda;
    //    double lambda=(mySpeciesArray(species)).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int sliceCounter=startSlice;sliceCounter<endSlice;sliceCounter+=skip){
    //      dVec r =mySpeciesArray(species,ptclNum,sliceCounter);
      //      dVec rp=mySpeciesArray(species,ptclNum,sliceCounter+skip);
      //      rpp=mySpeciesArray(species,ptclNum,sliceCounter+(skip>>1));
      dVec r = Path(sliceCounter,ptcl);
      dVec rp= Path(sliceCounter+skip,ptcl);
      dVec rpp=Path(sliceCounter+(skip>>1),ptcl);
      dVec rdiff = 
	DistanceTable->Velocity(sliceCounter, sliceCounter+skip, ptcl);
      ////      ///We've ignored boundary conditions here
      //dVec rbar=0.5*(r+rp);
      dVec rbar = r + 0.5*rdiff;
      dVec newDelta=GaussianRandomVec(sigma);

      for (int dim=0; dim<NDIM; dim++)
	{
	  while (newDelta[dim] > (0.5*Path.Box[dim]))
	    newDelta -= Path.Box[dim];
	  while (newDelta[dim] < (-(0.5*Path.Box[dim])))
	    newDelta += Path.Box[dim];
	}

      DistanceTable->PutInBox(newDelta);
      rpp=rbar+newDelta;
      logNewSampleProb=logNewSampleProb+
	(prefactorOfSampleProb-0.5*dot(newDelta,newDelta)/(sigma2));
      ////      ///Here we've stored the new position in the path
     //    mySpeciesArray.SetPos(species,ptclNum,sliceCounter+(skip>>1),rpp );
      Path.SetPos(sliceCounter+(skip>>1),ptcl,rpp);
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
  for (int ptcl=0;ptcl<particles.size();ptcl++){
    //    int species=particles(ptcl)(0);
    //    int ptclNum=particles(ptcl)(1);
    //    double lambda=mySpeciesArray(species).lambda;
    double lambda=Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0;//-NDIM/2.0*log(2*M_PI*sigma2);
    for (int sliceCounter=startSlice;sliceCounter<endSlice;sliceCounter+=skip){
      dVec r = Path(sliceCounter,ptcl);
      dVec rdiff = 
	DistanceTable->Velocity(sliceCounter, sliceCounter+skip, ptcl);
      dVec rp= Path(sliceCounter+skip,ptcl);
      dVec rpp=Path(sliceCounter+(skip>>1),ptcl);
      //      dVec r =mySpeciesArray(species,ptclNum,sliceCounter);
      //      dVec rp=mySpeciesArray(species,ptclNum,sliceCounter+skip);
      //      rpp    =mySpeciesArray(species,ptclNum,sliceCounter+(skip>>1));
      ///We've ignored boundary conditions here
      dVec rbar=r + 0.5*rdiff;
      dVec Delta= rpp - rbar;
      DistanceTable->PutInBox(Delta);
      logSampleProb=logSampleProb+
	(prefactorOfSampleProb-0.5*dot(Delta,Delta)/(sigma2));
    }
  }

  return logSampleProb;
}



#endif
