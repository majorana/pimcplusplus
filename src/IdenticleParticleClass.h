#ifndef IDENTICLE_PARTICLES_CLASS
#define IDENTICLE_PARTICLES_CLASS

#include "PathClass.h"

class PermutationClass
{



};


/*! This is an base class that holds all the information about
identical particles.  It may be specialized to hold specialized
information about particular types of particles. */
class IdenticalParticlesClass
{
public:
  int NumParticles;
  PathClass Path;
  PermutationClass Permutation;

  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classical particle.
  double lambda;
  
  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual double NodeAction (int Ptcl, int LinkNum) = 0;
};


class ElectronsClass : public IdenticalParticlesClass
{
public:
  double NodeAction (int Ptcl, int LinkNum);
  ElectronsClass();
  ~ElectronsClass();

};

#endif

