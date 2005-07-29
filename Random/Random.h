#ifndef RANDOM_CLASS
#define RANDOM_CLASS
#include <sprng.h>
#include "../MPI/Communication.h"


class RandomClass
{

private:
  int NumClones;
  int CloneNumber;
  int *LocalStream;
  int *CommonStream;
  CommunicatorClass &MyComm;
  int NumCommon, NumLocal;
public:
  ///Produces a double that is common among all processors
  inline double Common()
  { 
    double ranNum=sprng(CommonStream);  
    //    cerr<<"My random number(" << NumCommon << ") is " << ranNum 
    //	<< " and the stream is " << CommonStream 
    //	<<" and the proce num is " << MyComm.MyProc();
    //    cerr << "  NumLocal = " << NumLocal << endl;
    NumCommon++;
    return ranNum;
    //    return sprng(CommonStream);     
  }

  ///Produces a double that is unique to your processor
  inline double Local()
  {   
    NumLocal++;
    return sprng(LocalStream); 
  }

  ///Produces an int that is unique to your processor between 0 and max
  inline int LocalInt(int max)
  { return (int)floor(Local()*max); }
  
  /* normal random variate generator */
  // Returns a random number distributed according to
  // P(x) = 1/sqrt(2*pi * sigma) * exp ((x - mean)^2 /(2*sigma^2))
  inline double LocalGaussian(double sigma);
  
  /// Produces a guassian random vector with radius that has STD sigma
  inline void LocalGaussianVec(double sigma, Vec3 &c)
  { 
    c(0)=LocalGaussian(sigma); 
    c(1)=LocalGaussian(sigma); 
    c(2)=LocalGaussian(sigma); 
  }
  
  /// Produces a guassian random vector with radius that has STD sigma 
  inline void LocalGaussianVec(double sigma,Vec2 &c)
  { c(0)=LocalGaussian(sigma); c(1)=LocalGaussian(sigma); }
  
  inline double CommonGaussian(double sigma);
  /// Produces a guassian random vector with radius that has STD sigma
  inline void CommonGaussianVec(double sigma, Vec3 &c)
  { 
    c(0)=CommonGaussian(sigma); c(1)=CommonGaussian(sigma); c(2)=CommonGaussian(sigma); 
  }
  
  /// Produces a guassian random vector with radius that has STD sigma 
  inline void CommonGaussianVec(double sigma,Vec2 &c)
  { 
    c(0)=CommonGaussian(sigma); c(1)=CommonGaussian(sigma); 
  }


  ///Produces an int between 0 and max (not including max) that is
  ///common among processors  
  inline int CommonInt(int max)
  { 
    return (int)floor(Common()*max); 
  }
  
  
  void Init()
  {
    int seed = make_sprng_seed();
    //cerr<<"My seed is "<<seed<<endl;
    Init (seed);
  }
 
  void Init(int seed, int numClones=1)
  {
    cerr<<"I am "<<MyComm.MyProc()<<endl;
    cerr<<"I'm initing now"<<endl;
    NumClones=numClones;
    int myProc=MyComm.MyProc();
    int numProcs=MyComm.NumProcs();
    int procsPerClone=numProcs/NumClones;
    CloneNumber=myProc/procsPerClone;
    int commID=numProcs+CloneNumber;
    cerr<<"my data is "
	<<"my proc: "<<myProc<<" "
	<<"num clones: "<<NumClones<<" "
	<<"num procs: "<<numProcs<<" "
	<<"procs per clone: "<<procsPerClone<<" "
	<<"clone number: "<<CloneNumber<<" "
	<<"commID: "<<commID<<endl;
    ///LocalStream is a stream of random numbers unique to the node
    LocalStream = 
      init_sprng(SPRNG_DEFAULT, myProc,numProcs+NumClones,seed,SPRNG_DEFAULT);
    ///CommonStream is a stream of random numbers shared between all
    ///processors that are on the same clone 
    CommonStream = 
      init_sprng(SPRNG_DEFAULT, commID,numProcs+NumClones,seed,SPRNG_DEFAULT);

  }


  RandomClass(CommunicatorClass &comm) : 
    MyComm(comm), NumCommon(0), NumLocal(0)
  {

  }
};

inline double RandomClass::LocalGaussian(double sigma)
{                                      
  double V1, V2, S, X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft) {                   // Use second number from last computation
    X1 = X2;
    OneLeft = 0;
  }
  else {
    do {
      V1 = 2.0 * Local() - 1.0;
      V2 = 2.0 * Local() - 1.0;
      S = V1 * V1 + V2 * V2;
    } while ( S >= 1.0 );
    
    temp = sqrt( (-2.0 * log((double) S ) ) / S );
    X1 = V1 * temp;
    X2 = V2 * temp;
    OneLeft = 1;
  }  
  return( X1 * sigma );
}


inline double RandomClass::CommonGaussian(double sigma)
{                      
  double V1, V2, S, X1;
  static double X2;
  static int OneLeft = 0;
  double temp;

  if (OneLeft) {                   // Use second number from last computation
    X1 = X2;
    OneLeft = 0;
  }
  else {
    do {
      V1 = 2.0 * Common() - 1.0;
      V2 = 2.0 * Common() - 1.0;
      S = V1 * V1 + V2 * V2;
    } while ( S >= 1.0 );
    
    temp = sqrt( (-2.0 * log((double) S ) ) / S );
    X1 = V1 * temp;
    X2 = V2 * temp;
    OneLeft = 1;
  }  
  return( X1 * sigma );
}

#endif
