#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Common.h"
#include "SpeciesClass.h"
//#include "ActionClass.h"
#include "PathClass.h"
#include <Common/MPI/Communication.h>
#include "Actions/ActionsClass.h"
#include <Common/Random/Random.h>



/// This is the class that holds all of the information about the paths 
/// including all of the particles, the momized data, and the action.
/// Has routines for accepting and rejecting and for shifting data
/// between the processors.
class PathDataClass
{
private:
  int MyCloneNum;
  /// These store the time of start and the maximum wall time in seconds.
  int StartWallTime, MaxWallTime;
  int GetWallTime();
public:  
  int Seed;
  /// This defines a communicator for the group of processors working
  /// on this PathDataClass.
  /// This is for commmunication between nodes within a clone group.
  CommunicatorClass IntraComm;
  /// This is for communication between the rank 0 nodes of each clone
  /// group.  Hence, its between the clone groups.
  CommunicatorClass InterComm;
  /// This is the global MPIWORLD communicator.
  CommunicatorClass WorldComm;
  RandomClass Random;

  /// This function returns true if we have exceeded the maximum wall
  /// time.
  void SetMaxWallTime(int maxWallTime);
  bool ExceededWallTime();
  /// This object computes all actions.
  //  ActionClass Action; //(MemoizedDataClass,SpeciesArrayClass);
  ActionsClass Actions;
  /// The constructor that initializes the action and such
  int Join;
  PathClass Path;
  inline void ShiftData(int numTimeSlicesToShift){
    Path.ShiftData(numTimeSlicesToShift);
    Actions.ShiftData(numTimeSlicesToShift);
  }

  /// We are probaby going to have to move permutation
  /// information up here if we want it to notice
  /// the change in join.
  inline void MoveJoin(int newJoin)
  {
    Path.MoveJoin(Join,newJoin);
    Join=newJoin;
  }
  
  /// This function shifts the path to make the reference slice
  /// equal to the absolute slice position absSlice.  The join is left
  /// on the last slice, so the permutation is out of the way.
  void MoveRefSlice (int absSlice);

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
  inline SpeciesClass& Species(int species){
    return Path.Species(species);
  }


  inline int SpeciesNum(string speciesName);


  /// Returns the position of the particle of type species, particle
  /// number particle, and time slice timeSlice.
  inline dVec operator()(int timeSlice,int particle){
    return Path(timeSlice,particle);
  }

  /// Sets the position of the particle labeled by 
  /// (species, particle, timeSlice) to r and updates the time stamp
  /// of that piece of information.
  inline void SetPos(int timeSlice, int particle, const dVec& r){
    //    Path.SetPos(timeSlice,particle,r);
    Path(timeSlice,particle) = r;
  }

  inline int GetCloneNum() { return MyCloneNum; }

  void Read (IOSectionClass &in);
  PathDataClass() : 
    /* Action(*this), */Actions(*this), Random(WorldComm), 
			Path(IntraComm,Random, Actions), MaxWallTime(-1)
  { 
    Join = 1; 
    StartWallTime = GetWallTime();
  }


};
 
inline int PathDataClass::SpeciesNum(string name)
{

  for (int spec=0;spec<NumSpecies();spec++){
    if (Species(spec).Name==name){
      return spec;
    }
  }
  return -1;
}

inline void PathDataClass::AcceptMove(int startTimeSlice,int endTimeSlice,
			       const Array <int,1> &activeParticles)
{
  Actions.AcceptCopy (startTimeSlice, endTimeSlice, activeParticles);
  Path.AcceptCopy(startTimeSlice,endTimeSlice,activeParticles);
}

inline void PathDataClass::RejectMove(int startTimeSlice,int endTimeSlice,
			       const Array <int,1> &activeParticles)
{
  Actions.RejectCopy (startTimeSlice, endTimeSlice, activeParticles);
  Path.RejectCopy(startTimeSlice,endTimeSlice,activeParticles); 
}


#endif
