#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Common.h"
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
  CommClass Communicator;
  int NumTimeSlices;
  void acceptMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice);
  void rejectMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice);
};
#endif
