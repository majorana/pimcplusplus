#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Common.h"
#include "IdenticalParticlesClass.h"
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
  CommunicatorClass Communicator;
  int NumTimeSlices;
  void acceptMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice);
  void rejectMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice);
  inline int NumSpecies(){
    return IdenticalParticleArray.size();
  }
  inline PathClass& operator()(int Species){
    return IdenticalParticleArray(Species).Path;
  }

  inline dVec operator()(int Species,int Particle,int TimeSlice){
    return IdenticalParticleArray(Species,Particle,TimeSlice);
  }
  inline void SetPos(int Species, int Particle, int TimeSlice,const dVec& r){
    IdenticalParticleArray.SetPos(Species,Particle,TimeSlice,r);
  }

};

#endif
