#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "CubicSpline.h"
#include "IdenticleParticleClass.h"
#include "MemoizedDataClass.h"

class SavedPairActionClass
{



};

/*! $[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n \sum_{j=1}^k
u_{kj}(q;\tau)z^{2j}s^{2(k-j)} */

class PairActionClass
{

private:
  Array<double,1> tempukjArray;

public:
  /// This stores the coefficients of the expansion specified above.
  /// The array index is over the levels.  You call it with the q
  /// value and a temporary array to get all of the values in that
  /// column. 
  Array<MultiCubicSpline,1> ukj; //(level )
  inline double calcUrrptau(double s,double q,double z,int level);
  int n;
  double tau;

  
  
  
  
};

inline double PairActionClass::calcUrrptau(double s,double q,double z,int level)
{
  double sum=0;
  double r=0.5*(q+z);
  double rprime=0.5*(q-z);
///I'm about to change this line to make it work  sum=sum+(ukj(level,0))(r)+(ukj(level,0))(rprime);//this is the endpoint action
  sum=sum+(ukj(level))(0,r)+(ukj(level))(0,rprime);//this is the endpoint action
  double zsquared=z*z;
  double ssquared=s*s;
  double ssquaredinverse=1/ssquared;
  double Sto2k=2;
  (ukj(level))(q,tempukjArray); 
  for (int k=1;k<=n;k++){  
				
    double Zto2j=1;
    double currS=Sto2k;

    for (int j=0;j<=k;j++){

      double cof=tempukjArray(k*(k+1)/2+j); //indexing into the 2darray
      sum=sum+cof*Zto2j*currS;
      
      
      Zto2j*=zsquared;
      currS=currS*ssquaredinverse;				
    }				
    Sto2k=Sto2k*ssquared;
  }
  
  return sum; //I hope this is the right thing to return 
}

/*! This is the class that controls all of the actions and is in
  charge of calculating them. When this is initialized a pointer needs
  to be sent that has the memoizedData and IdenticalParticleClass */ 

class ActionClass
{


private:
public:
  Array<PairActionClass,1> PairActionVector;
  Array<int,2> PairMatrix;
  Array<SavedPairActionClass,2> SavedPairActionArray;
  Array<IdenticalParticlesClass,1> *myIdenticalParticleArray;
  double tau;
  bool calcTotalAction(Array<ParticleID,1> changedParticles,int level);
  MemoizedDataClass *myMemoizedData;
  inline void SampleParticles(Array<ParticleID,1> particles,int startSlice,int endSlice,int level,double&, double&);
  void calcTotalAction();



};

/*! 
Samples the particles in "particles" and the time slices between
startSlice and endSlice (not inclusive). Samples from a free gaussian
of variance \$\sigma=\sqrt{2*\lambda*\tau} \$. Return 
the logarithm of the total sampling probability of the density (i.e \$
\frac{-NDIM}{2}*\log{2*\pi*\sigma^2}-0.5*\frac{\delta \cdot
\delta}{\sigma * \sigma}  \$ )

!*/

inline void ActionClass::SampleParticles(Array<ParticleID,1> particles,int startSlice,int endSlice, int level,double &logNewSampleProb, double &logOldSampleProb)
{
  dVec rpp;
  int skip = 1<<(level+1);
  logNewSampleProb=1;
  logOldSampleProb=1;

  double levelTau = 0.5*tau*skip;
  for (int ptcl=0;ptcl<particles.size();ptcl++){
    int species=particles(ptcl)(0);
    int ptclNum=particles(ptcl)(1);
    double lambda=((*myIdenticalParticleArray)(species)).lambda;
    double sigma2=(2*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=-NDIM/2*log(2*M_PI*sigma2);
    for (int sliceCounter=startSlice;sliceCounter<endSlice;sliceCounter+=skip){
      dVec r =(*myIdenticalParticleArray)(species).Path(ptclNum,sliceCounter);
      dVec rp=(*myIdenticalParticleArray)(species).Path(ptclNum,sliceCounter+skip);
      rpp=(*myIdenticalParticleArray)(species).Path(ptclNum,sliceCounter+skip>>1);
      ///We've ignored boundary conditions here
      dVec rbar=0.5*(r+rp);
      dVec newDelta=GuassianRandomVec(sigma);
      
      dVec oldDelta=rpp-rbar;
      rpp=rbar+newDelta;
      logNewSampleProb=logNewSampleProb+(prefactorOfSampleProb-0.5*dot(newDelta,newDelta)/(sigma2));
      logOldSampleProb=logOldSampleProb+(prefactorOfSampleProb-0.5*dot(oldDelta,oldDelta)/(sigma2));
      ///Here we've stored the new position in the path
      (*myIdenticalParticleArray)(species).Path(ptclNum,sliceCounter+skip>>1).setPos(rpp);
      
    }
  }


}

#endif
