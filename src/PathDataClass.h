#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "IdenticleParticleClass.h"
#include "MemoizedDataClass.h"
#include "ActionClass.h"

/*! This is the class that holds all of the information about the paths including all of the particles, the momized data, and the action */

class PathDataClass
{
public:
  ActionClass TotalAction;
  MemoizedDataClass MemoizedData;
  ArrayOfIdenticalParticlesClass IdenticalParticleArray;
    //  Array<IdenticalParticleClass,1> IdenticalParticleArray;

};
#endif
