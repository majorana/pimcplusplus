#ifndef PATH_DATA_CLASS
#define PATH_DATA_CLASS

#include "Common.h"
#include "SpeciesClass.h"
#include "ActionClass.h"
#include "PathClass.h"
#include "CommunicatorClass.h"
#include "DistanceTableClass.h"
#include "BCClass.h"

#include "Common/Random/Random.h"
//#include "DistanceTableFreeClass.h"

/// This is the class that holds all of the information about the paths 
/// including all of the particles, the momized data, and the action.
/// Has routines for accepting and rejecting and for shifting data
/// between the processors.
class PathDataClass
{
public:
  BCClass BC;
  DistanceTableClass *DistanceTable;
  

  /// This defines a communicator for the group of processors working
  /// on this PathDataClass.
  PIMCCommunicatorClass Communicator;

  /// This object computes all actions.
  ActionClass Action; //(MemoizedDataClass,SpeciesArrayClass);

  /// The constructor that initializes the action and such
  int Join;
  PathClass Path;
  inline void ShiftData(int numTimeSlicesToShift){
    Path.ShiftData(numTimeSlicesToShift);
    DistanceTable->ShiftData(numTimeSlicesToShift, Communicator);
  }

  /// We are probaby going to have to move permutation
  /// information up here if we want it to notice
  /// the change in join.
  inline void MoveJoin(int newJoin)
  {
    Path.MoveJoin(Join,newJoin);
    DistanceTable->MoveJoin(Join,newJoin);
    Join=newJoin;
  }
  
  inline void Update(int slice,const Array<int,1> &ptclArray)
  {
    DistanceTable->Update(slice,ptclArray);
  }

  inline void Update(int startSlice,int endSlice,
		     const Array <int,1> &ptclArray,int level)
  {
    int skip = 1<<(level+1);
    int skipo2 = skip >> 1;
    for (int slice=startSlice;slice<endSlice;slice+=skip){
      Update(slice+skipo2,ptclArray);
    }
  }

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
    Path.SetPos(timeSlice,particle,r);
  }

  inline PathDataClass() : 
    Action(Path), Path(Communicator), BC(Path)
  { Join = 1; }

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
  Path.AcceptCopy(startTimeSlice,endTimeSlice,activeParticles);
  DistanceTable->AcceptCopy(startTimeSlice,endTimeSlice,activeParticles);
}

inline void PathDataClass::RejectMove(int startTimeSlice,int endTimeSlice,
			       const Array <int,1> &activeParticles)
{
  Path.RejectCopy(startTimeSlice,endTimeSlice,activeParticles);  
  DistanceTable->RejectCopy(startTimeSlice,endTimeSlice,activeParticles);
}


#endif
