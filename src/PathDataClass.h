#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Common.h"
#include "IdenticalParticlesClass.h"
#include "MemoizedDataClass.h"
#include "ActionClass.h"

/// This is the class that holds all of the information about the paths 
/// including all of the particles, the momized data, and the action.
/// Has routines for accepting and rejecting and for shifting data
/// between the processors.
class PathDataClass
{
public:
  /// This object computes all actions.
  ActionClass TotalAction;
  /// This object holds computed data which is used multiple times.
  MemoizedDataClass MemoizedData;
  /// This object functions as an array of SpeciesClass objects.
  ArrayOfIdenticalParticlesClass IdenticalParticleArray;
  /// This defines a communicator for the group of processors working
  /// on this PathDataClass.
  CommunicatorClass Communicator;
  /// Returns the number of time slices.
  inline int NumTimeSlices()
  {  return (IdenticalParticleArray(0).NumTimeSlices());  }
  /// Do all copies necessary to accept a move.
  void AcceptMove(Array <ParticleID,1> activeParticles,
		  int startTimeSlice,int endTimeSlice);
  /// Do all copies necessary to accept a move.
  void RejectMove(Array <ParticleID,1> activeParticles,
		  int startTimeSlice,int endTimeSlice);
  /// Returns the number of particle species in the path
  inline int NumSpecies(){
    return IdenticalParticleArray.Size();
  }
  /// Returns a reference to the SpeciesClass object of number species
  inline IdenticalParticlesClass& operator()(int species){
    return IdenticalParticleArray(species);
  }
  /// Returns the position of the particle of type species, particle
  /// number particle, and time slice timeSlice.
  inline dVec operator()(int species,int particle,int timeSlice){
    return IdenticalParticleArray(species,particle,timeSlice);
  }
  /// Sets the position of the particle labeled by 
  /// (species, particle, timeSlice) to r and updates the time stamp
  /// of that piece of information.
  inline void SetPos(int species, int particle, int timeSlice,const dVec& r){
    IdenticalParticleArray.SetPos(species,particle,timeSlice,r);
  }
};

#endif
