#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "CubicSpline.h"


class PairActionClass
{
  


}

/*! This is the class that controls all of the actions and is in
  charge of calculating them. When this is initialized a pointer needs
  to be sent that has the memoizedData and IdenticalParticleClass */ 

class ActionClass
{
public:
  Array<PairActionClass,1> PairActionVector;
  Array<int,2> PairMatrix;
  Array<SavedPairActionClass,2> SavedPairActionArray;
  IdenticalParticleClass *myIdenticalParticleClass;
  MemoizedDataClass *myMemoizedDataClass;
  calcTotalAction();
private:


};

