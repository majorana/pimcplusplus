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
  /// Stores the permutation for my set of time-slices
  PermutationClass Permutation;
  
  inline int NumParticles()
  { return Path.NumParticles(); }

  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classical particle.
  double lambda;
  
  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual double NodeAction (int Ptcl, int LinkNum) = 0;
};


class ElectronsClass : public SpeciesClass
{
public:
  double NodeAction (int Ptcl, int LinkNum);
  ElectronsClass();
  ~ElectronsClass();

};

class ProtonsClass : public SpeciesClass
{
public:
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

