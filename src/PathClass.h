#ifndef PATH_CLASS_H
#define PATH_CLASS_H

#include "Common/IO/InputOutput.h"
#include "MirroredArrayClass.h"
#include "SpeciesClass.h"
#include "Common/Random/Random.h"

///The number of time slices is the number of slices on this processor.
///In all cases this processor shares a time slice with the processor 
///ahead of it and behind it. The convention for the shared slices
///is that the processor owns its first but not its last slice.
class PathClass
{
private:
  /// Path stores the position of all the particles at all time
  /// slices.  The order for access is timeslice, particle
  MirroredArrayClass<dVec> Path;
  /// Stores what species a particle belongs to
  Array<int,1> SpeciesNumber;
  Array<SpeciesClass *,1> SpeciesArray;
  int MyNumSlices;

  PIMCCommunicatorClass &Communicator;

public:
  // True if we're using periodic boundary conditions.
  bool UsePBC;
  RandomClass Random;
  MirroredArrayClass1D<int> Permutation;
  inline void Print(){Path.Print();}
  void Read(IOSectionClass &inSection);
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

  inline int SpeciesNum (string name)
  {
    int i=0;
    while ((i<SpeciesArray.size()) && (SpeciesArray(i)->Name != name))
      i++;
    if (i == SpeciesArray.size())
      return -1;
    else
      return i;
  }


  /// Shifts the data to other processors or to yourself if there 
  /// are no other processors
  inline void ShiftData(int sliceToShift)
  {Path.ShiftData(sliceToShift,Communicator);}

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
  int TotalNumSlices;
  /// Returns the number of particle Species
  inline int NumSpecies() {return SpeciesArray.size();}
  inline int NumParticles() { return Path.NumParticles();}
  ///The number of time slices is the number of slices on this processor.
  ///In all cases this processor shares a time slice with the processor 
  ///ahead of it and behind it. The convention for the shared slices
  ///is that the processor owns its first but not its last slice.
  inline int NumTimeSlices() { return Path.NumTimeSlices();}

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
    assert(TotalNumSlices>0);
    int myProc=Communicator.MyProc();
    int numProcs=Communicator.NumProcs();
    ///Everybody gets the same number of time slices if possible.
    ///Otherwise the earlier processors get the extra one slice 
    ///until we run out of extra slices.
    ///The last slice on processor i is the first slices on processor i+1
    MyNumSlices=TotalNumSlices/numProcs+1+(myProc<(TotalNumSlices % numProcs));
    cerr<<"Numprocs is "<<numProcs<<endl;
    cerr<<"mynumslices: "<<MyNumSlices<<endl;
    

    int numParticles = 0;
    /// Set the particle range for the new species
    for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
      SpeciesArray(speciesNum)->FirstPtcl = numParticles;
      numParticles=numParticles + SpeciesArray(speciesNum)->NumParticles;
      SpeciesArray(speciesNum)->LastPtcl= numParticles-1;
    }
    Path.Resize(MyNumSlices,numParticles);
    cerr<<"I've resized to "<<MyNumSlices<<" "<<numParticles<<endl;
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

  PathClass(PIMCCommunicatorClass &communicator): Communicator(communicator), Random(Communicator)
    {
      //      NumSpecies = 0;
      TotalNumSlices=0;
      Random.Init();
    }
};

#endif
