#ifndef SPECIES_CLASS
#define SPECIES_CLASS

#include "PathClass.h"
#include "PermutationClass.h"

/// This is an base class that holds all the information about
/// identical particles.  It may be specialized to hold specialized
/// information about particular types of particles.
class SpeciesClass
{
public:
  /// Stores the positions and timestamp for all particles and TimeSlices
  PathClass Path;
  int Join;
  /// Stores the permutation for my set of time-slices. This needs to be resized at some point!!
  MirroredArrayClass1D<int> Permutation;
  
  
  inline int NumParticles()
  { return Path.NumParticles(); }
  inline int NumTimeSlices()
  { return Path.NumTimeSlices(); }


  ////This function moves the join from one place to another.
  ///Do not move the join to slice 0 or slice n. I think this 
  ///will cause problems. 

  inline void MoveJoin(int newJoin)
  {
    Path.MoveJoin(Permutation,Join,newJoin);
    Join=newJoin;
    return;
  }


  void ShiftData(int sliceToShift, CommunicatorClass &communicator);
  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classical particle.
  double lambda;
  
  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual double NodeAction (int Ptcl, int LinkNum) = 0;
};




///This is the inherited class that holds the information about
///the electrons. Eventually this will probably be turned into
///a fermion class.  Currently, no actual information about the electrons
///are actually included.
class ElectronsClass : public SpeciesClass
{
public:
  ///When we make this work, this calculates the NodeActions
  double NodeAction (int Ptcl, int LinkNum);
  ElectronsClass();
  ~ElectronsClass();

};


///This is the inherited class that holds the information about
///the protons.  Currently no useful information about protons is contained
///here
class ProtonsClass : public SpeciesClass
{
public:
  ///Just returns 0 until we do something more intelligble.
  double NodeAction (int Ptcl, int LinkNum)
    {
      return (0.0);
    }  
  ProtonsClass()
    {
    }
  ~ProtonsClass()
    {
    }
  
};

#endif

