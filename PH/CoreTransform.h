#ifndef CORE_TRANSFORM_H
#define CORE_TRANSFORM_H

#include "PH.h"

/// The core transform class describes a coordinate transformation
/// used in calculating the pair density matrix through matrix
/// squaring for pseudohamiltonians .  It maps between real radial
/// distances, r, and transformed distances, x.  This transformation
/// allows the radial equation for the pseudohamiltonians to have the
/// same canonical form as a local radial potential, thereby allowing
/// us to use the traditional matrix squaring formalism.
class CoreTransform
{
  scalar rMax, xMax;
  LinearGrid xgrid, rgrid;
  CubicSpline x_of_r;
  CubicSpline r_of_x;

public:
  /// Convert from x coordinate to r coordinate
  inline scalar x2r(scalar x)
  {
    if (x >= xMax)
      return ((x-xMax)+rMax);
    else
      return (r_of_x(x));
  }

  /// Convert from r coordinate to x coordinate
  inline scalar r2x(scalar r)
  {
    if (r >= rMax)
      return ((r-rMax)+xMax);
    else
      return (x_of_r(r));
  }
  
  /// Calculate the transform from the pseudohamiltonian.
  void Initialize(PseudoHamiltonian *PH, int NumPoints);
};

#endif
