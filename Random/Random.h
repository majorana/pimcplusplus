#ifndef RANDOM_CLASS
#define RANDOM_CLASS
#include "sprng.h"
#include "../MPI/Communication.h"


class RandomClass
{

private:
  int NumClones;
  int CloneNumber;
  int *LocalStream;
  int *CommonStream;
  CommunicatorClass &MyComm;

public:
  ///Produces a double that is common among all processors
  inline double Common()
  {   return sprng(CommonStream);     }

  ///Produces a double that is unique to your processor
  inline double Local()
  {   return sprng(LocalStream); }
  
  /* normal random variate generator */
  // Returns a random number distributed according to
  // P(x) = 1/sqrt(2*pi * sigma) * exp ((x - mean)^2 /(2*sigma^2))
  inline double LocalGaussian(double mean, double sigma);                                      
    
  /// Produces a guassian random vector with radius that has variance
  /// sigma 
  inline void LocalGaussianVec(double sigma, Vec3 &c)
  {
    c(0)=LocalGaussian(0.0,sigma);
    c(1)=LocalGaussian(0.0,sigma);
    c(2)=LocalGaussian(0.0,sigma);    
  }
  
  /// Produces a guassian random vector with radius that has variance
  /// sigma 
  inline void LocalGaussianVec(double sigma,Vec2 &c)
  {
    c(0)=LocalGaussian(0.0,sigma);
    c(1)=LocalGaussian(0.0,sigma);
  }
  
  
  void Init()
  {
    int seed = make_sprng_seed();
    cerr<<"My seed is "<<seed<<endl;
    Init (seed);
  }

  void Init(int seed)
  {
    int myProc=MyComm.MyProc();
    int numProcs=MyComm.NumProcs();
    int procsPerClone=numProcs/NumClones;
    CloneNumber=myProc/procsPerClone;
    int commID=numProcs+CloneNumber;
    LocalStream = init_sprng(SPRNG_DEFAULT, myProc,numProcs+NumClones,seed,SPRNG_DEFAULT);
    CommonStream = init_sprng(SPRNG_DEFAULT, commID,numProcs+NumClones,seed,SPRNG_DEFAULT);    
  }


  RandomClass(CommunicatorClass &comm, int numClones=1) : MyComm(comm)
  {
    NumClones=numClones;
  }
};

inline double RandomClass::LocalGaussian(double mean, double sigma)
{                                      
  double V1, V2, S, X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft)                   // Use second number from last computation
    {
      X1 = X2;
      OneLeft = 0;
    }
  else
    {
      do
        {
          V1 = 2.0 * Local() - 1.0;
          V2 = 2.0 * Local() - 1.0;
          S = V1 * V1 + V2 * V2;
        } while ( S >= 1.0 );
      
      temp = sqrt( (-2.0 * log((double) S ) ) / S );
      X1 = V1 * temp;
      X2 = V2 * temp;
      OneLeft = 1;
    }  
  return( mean + X1 * sigma );
}

#endif
