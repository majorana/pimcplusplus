#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "IdenticalParticleClass.h"


/*! This is the class that holds all of the information about the paths including all of the particles, the momized data, and the action */

class PAthDataClass
{
public:
  ActionClass TotalAction;
  MemoizedDataClass MemoizedData;
  Array<IdenticalParticleClass,1> IdenticalParticleArray;

};
