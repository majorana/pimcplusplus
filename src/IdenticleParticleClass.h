#ifndef IDENTICLE_PARTICLES_CLASS
#define IDENTICLE_PARTICLES_CLASS

#include "PathClass.h"




/*! This is an base class that holds all the information about
identicle particles.  It may be specialized to hold specialized
information about particular types of particles. */
class IdenticleParticlesClass
{
public:
  PathClass Path;
  PermutationClass Permutation;

  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classicle particle.
  double lambda;
  
  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual double NodeAction (int Ptcl, int LinkNum);
}


#endif
