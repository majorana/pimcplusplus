#include "Common.h"





int GetCurrentTimeStamp()
{
  return 0;
}



void setMode(ModeType theMode)
{

  if (theMode==MOVEMODE){
    Write1=0;
    Write2=0;
  }

  else  if (theMode==OBSERVABLEMODE){
    Write1=0;
    Write2=1;
  }


}

dVec GuassianRandomVec(double sigma)
{
  dVec c;
  c(0)=0;
  c(1)=1;
  c(2)=2;
  return c;

}

dVec dVecSubtract(dVec c,dVec d)
{

  dVec b;
  b(0)=0;
  b(1)=1;
  b(2)=2;
  return b;

}
double distSqrd(dVec a, dVec b) //Did I do this right?
{
  dVec Dist= distSqrd(a,b);
  return (Dist(0)*Dist(0)+Dist(1)*Dist(1)+Dist(2)*Dist(2));
}

  
