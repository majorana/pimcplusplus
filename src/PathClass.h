#ifndef PATH_CLASS_H
#define PATH_CLASS_H


#include "MirroredArrayClass.h"
#include "SpeciesClass.h"

class PathClass
{
private:
  /// Path stores the position of all the particles at all time
  /// slices.  The order for access is timeslice, particle
  MirroredArrayClass<dVec> Path;
  Array<int,1> SpeciesNumber;
  Array<SpeciesClass *,1> SpeciesArray;
  int TimeSliceNumber;
public:
  MirroredArrayClass1D<int> Permutation;
  ///A scratch array to hold a boolean indicating whether we've
  ///looped over this particle yet
  Array<bool,1> DoPtcl;
  dVec Box;
  inline void MoveJoin(int oldJoin, int newJoin){
    Path.MoveJoin(Permutation,oldJoin,newJoin);
  }
      
  inline void AcceptCopy(int startTimeSlice,int endTimeSlice, 
	     const Array <int,1> &activeParticle)
  {
      Path.AcceptCopy(startTimeSlice,endTimeSlice,activeParticle);

  }

  inline void RejectCopy(int startTimeSlice,int endTimeSlice, 
	     const Array <int,1> &activeParticle )
  {
      Path.RejectCopy(startTimeSlice,endTimeSlice,activeParticle);

  }




  /// Shifts the data to other processors or to yourself if there 
  /// are no other processors
  inline void ShiftData(int sliceToShift, CommunicatorClass &communicator)
  {Path.ShiftData(sliceToShift,communicator);}

  /// Return what species type a particle belongs to;
  inline int ParticleSpeciesNum(int ptcl)
  {
    return (SpeciesNumber(ptcl));
  }
  /// Return a reference to the species that particle ptcl belongs to
  inline SpeciesClass& ParticleSpecies(int ptcl) 
  {
    return *(SpeciesArray(SpeciesNumber(ptcl)));
  }
  /// Return a species references
  inline SpeciesClass& Species(int speciesNum)
  {
    return (*(SpeciesArray(speciesNum)));
  }
  /// Returns the number of particle Species
  inline int NumSpecies() {return SpeciesArray.size();}
  inline int NumParticles() { return Path.NumParticles();}
  inline int NumTimeSlices() { return Path.NumTimeSlices();}
  inline void SetTimeSlices(int tSlices){TimeSliceNumber=tSlices;}
  /// Returns the position of particle ptcl at time slice timeSlice
  inline dVec operator() (int timeSlice, int ptcl)
  { return Path(timeSlice, ptcl); }
  /// Set the position of particle ptcl at time slice timeSlice
  inline void SetPos (int timeSlice, int ptcl, dVec r)
  { Path.Set(timeSlice, ptcl, r); }

  void AddSpecies (SpeciesClass *newSpecies)
  {
    int numSpecies = SpeciesArray.size();
    /// Add an element for the new species
    SpeciesArray.resizeAndPreserve(numSpecies+1);
    SpeciesArray(numSpecies) = newSpecies;
  }

  void Allocate()
  {
    assert(TimeSliceNumber>0);
    int numParticles = 0;
    /// Set the particle range for the new species
    for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
      SpeciesArray(speciesNum)->FirstPtcl = numParticles;
      numParticles=numParticles + SpeciesArray(speciesNum)->NumParticles;
      SpeciesArray(speciesNum)->LastPtcl= numParticles-1;
    }
    Path.Resize(TimeSliceNumber,numParticles);
    Permutation.Resize(numParticles);
    SpeciesNumber.resize(numParticles);
    DoPtcl.resize(numParticles);
    /// Assign the species number to the SpeciesNumber array
    for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
      for (int i=SpeciesArray(speciesNum)->FirstPtcl; 
	   i<= SpeciesArray(speciesNum)->LastPtcl; i++)
	SpeciesNumber(i) = speciesNum;
    }
    //Sets to the identity permutaiton 
    for (int ptcl=0;ptcl<Permutation.NumParticles();ptcl++){
      Permutation.Set(ptcl,ptcl);
    }
  }

  PathClass()
    {
      //      NumSpecies = 0;

      TimeSliceNumber=0;
      
    }
};

#endif
