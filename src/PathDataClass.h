#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Common.h"
#include "SpeciesClass.h"
#include "MemoizedDataClass.h"
#include "ActionClass.h"

/// This is the class that holds all of the information about the paths 
/// including all of the particles, the momized data, and the action.
/// Has routines for accepting and rejecting and for shifting data
/// between the processors.
class PathDataClass
{
private:
  int Join;

public:
  ///The constructor that initializes the action and such

  PathClass Path;
  inline PathDataClass() : Action(MemoizedData,SpeciesArray)
    {

      Join=1;
  
    }

  DistanceTableClass DistanceTable;
  /// This object functions as an array of SpeciesClass objects.

  /// This defines a communicator for the group of processors working
  /// on this PathDataClass.
  CommunicatorClass Communicator;

  /// This object computes all actions.
  ActionClass Action; //(MemoizedDataClass,SpeciesArrayClass);


  /// Returns the number of time slices.
  inline int NumTimeSlices()
  {  return Path.NumTimeSlices();  }
  /// Do all copies necessary to accept a move.
  inline void AcceptMove(int startTimeSlice,int endTimeSlice,
			 const Array <int,1> &activeParticles);
		  
  /// Do all copies necessary to accept a move.
  inline void RejectMove(int startTimeSlice,int endTimeSlice,
			 const Array <int,1> &activeParticles);
		  
  /// Returns the number of particle species in the path
  inline int NumSpecies(){
    return Path.NumSpecies();
  }
  /// Returns the total number of particles
  inline int NumParticles(){
    return Path.NumParticles();
  }
  /// Returns a reference to the SpeciesClass object of number species
  //  inline SpeciesClass& operator()(int species){
  //    return Path.SpeciesArray(species);
  //  }
  /// Returns the position of the particle of type species, particle
  /// number particle, and time slice timeSlice.
  inline dVec operator()(int timeSlice,int particle){
    return Path(timeSlice,particle);
  }
  /// Sets the position of the particle labeled by 
  /// (species, particle, timeSlice) to r and updates the time stamp
  /// of that piece of information.
  inline void SetPos(int timeSlice, int particle, const dVec& r){
    SpeciesArray.SetPos(timeSlice,particle,r);
  }
};


inline void PathDataClass::AcceptMove(int startTimeSlice,int endTimeSlice,
			       const Array <int,1> &activeParticles)
{

  Path.AcceptCopy(startTimeSlice,endTimeSlice,activeParticles);
  DistanceTable.AcceptCopy(startTimeSlice,endTimeSlice,activeParticles);
  
}

inline void PathDataClass::RejectMove(int startTimeSlice,int endTimeSlice,
			       const Array <ParticleID,1> &activeParticles)
{
  Path.RejectCopy(startTimeSlice,endTimeSlice,activeParticles);  
  DistanceTable.RejectCopy(startTimeSlice,endTimeSlice,activeParticles);

  
}


#endif
