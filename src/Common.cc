#include "Common.h"


/* normal random variate generator */
// Returns a random number distributed according to
// P(x) = 1/sqrt(2*pi * sigma) * exp ((x - mean)^2 /(2*sigma^2))
scalar GaussianRandNum(scalar mean, scalar sigma)
{                                      
  scalar V1, V2, S, X1;
  static scalar X2;
  static int OneLeft = 0;
  scalar temp;

  if (OneLeft)                   // Use second number from last computation
    {
      X1 = X2;
      OneLeft = 0;
    }
  else
    {
      do
        {
          V1 = 2.0 * sprng() - 1.0;
          V2 = 2.0 * sprng() - 1.0;
          S = V1 * V1 + V2 * V2;
        } while ( S >= 1.0 );
      
      temp = sqrt( (-2.0 * log((double) S ) ) / S );
      X1 = V1 * temp;
      X2 = V2 * temp;
      OneLeft = 1;
    }  
  return( mean + X1 * sigma );
}






int GetCurrentTimeStamp()
{
  return 0;
}


/// In NEW mode, write only to the first copy,in OLD MODE 
/// write to the section copy,  in OBSERVABLE
/// mode, write to both copies.

void setMode(ModeType theMode)
{

  if (theMode==OLDMODE){
    Write1=1;
    Write2=1;
  }
  else if (theMode==NEWMODE)
    {
      Write1=0;
      Write2=0;
    }
  else  if (theMode==BOTHMODE){
    Write1=0;
    Write2=1;
  }


}

dVec GaussianRandomVec(double sigma)
{
  dVec c;
  c(0)=GaussianRandNum(0.0,sigma);
  c(1)=GaussianRandNum(0.0,sigma);
  c(2)=GaussianRandNum(0.0,sigma);
  return c;

}




double distSqrd(dVec a, dVec b) //Did I do this right?
{
  dVec Disp = a - b;
  return (dot(Disp,Disp));
}

  
