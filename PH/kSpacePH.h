#ifndef K_SPACE_PH_H
#define K_SPACE_PH_H

#include "../Blitz.h"

#include "PotentialBase.h"

/// kSpacePH class takes a PH as its constructor argument.  Its only
/// public member function, V, returns the k-space component.  The
/// folumae used in this class come from 
/// Foulkes and Schluter Phys. Rev. B 42, 11505 (1990)
class kSpacePH
{
protected:
  Potential &PH;
  double Vk (double k);
  double a  (double k);
  /// The perpendicular component of the b tensor
  double bPerp (double k);
  /// The perpendicular component of the b tensor
  double bPar  (double k);
  /// Returns the F tensor for k
  TinyMatrix<double,3,3> F (double k);
public:
  double V (Vec3 k, Vec3 G, Vec3 Gp);

  kSpacePH (Potential &ph) : PH(ph)
  {
    // do nothing else for now
  }
};

#endif
