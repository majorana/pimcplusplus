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
  /// The tail is fit to the potential to the form: 
  /// V = Ctail/r + Ctail2/r^2 + Ctail3/r^3
  double Ctail1, Ctail2, Ctail3;

  /// Member functions 
  double Vk (double k);
  double a  (double k);
  /// The perpendicular component of the b tensor
  double bPerp (double k);
  /// The perpendicular component of the b tensor
  double bPar  (double k);
  /// Returns the F tensor for k
  TinyMatrix<double,3,3> Ftensor (Vec3 deltaG);
  bool HaveTailCoefs;
  double R1, R2;
public:
  /// Calculates the values of Ctail1-3.  r1 and r2 specify the start
  /// and end of the region of the potential fit to the form given
  /// above. 
  void CalcTailCoefs(double r1, double r2);

  /// This returns the k-space representation of the
  /// pseudohamiltonian.  Note that it does not contain the 1/Vcell
  /// volume factor, and that you MUST call CalcTailCoefs before
  /// calling this fuction.
  double V (Vec3 k, Vec3 G, Vec3 Gp);
  

  kSpacePH (Potential &ph) : PH(ph), HaveTailCoefs(false)
  {
    // do nothing else for now
  }
};

#endif
