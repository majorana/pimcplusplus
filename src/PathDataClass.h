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
  ActionClass Action;
  /// This object holds computed data which is used multiple times.
  MemoizedDataClass MemoizedData;
  /// This object functions as an array of SpeciesClass objects.
  SpeciesArrayClass SpeciesArray;
  /// This defines a communicator for the group of processors working
  /// on this PathDataClass.
  CommunicatorClass Communicator;
  /// Returns the number of time slices.
  inline int NumTimeSlices()
  {  return (SpeciesArray(0).NumTimeSlices());  }
  /// Do all copies necessary to accept a move.
  void AcceptMove(Array <ParticleID,1> activeParticles,
		  int startTimeSlice,int endTimeSlice);
  /// Do all copies necessary to accept a move.
  void RejectMove(Array <ParticleID,1> activeParticles,
		  int startTimeSlice,int endTimeSlice);
  /// Returns the number of particle species in the path
  inline int NumSpecies(){
    return SpeciesArray.Size();
  }
  /// Returns a reference to the SpeciesClass object of number species
  inline IdenticalParticlesClass& operator()(int species){
    return SpeciesArray(species);
  }
  /// Returns the position of the particle of type species, particle
  /// number particle, and time slice timeSlice.
  inline dVec operator()(int species,int particle,int timeSlice){
    return SpeciesArray(species,particle,timeSlice);
  }
  /// Sets the position of the particle labeled by 
  /// (species, particle, timeSlice) to r and updates the time stamp
  /// of that piece of information.
  inline void SetPos(int species, int particle, int timeSlice,const dVec& r){
    SpeciesArray.SetPos(species,particle,timeSlice,r);
  }
};

#endif
