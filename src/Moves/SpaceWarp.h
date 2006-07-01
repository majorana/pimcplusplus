#ifndef SPACE_WARP_H
#define SPACE_WARP_H

#include <Common/Blitz.h>

class SpaceWarpClass
{
protected:
  int N;
  /// The current ion positions
  Array<Vec3,1> Rions;
  /// The proposed change in ion positions
  Array<Vec3,1> DeltaRions;
  Array<Vec3,1> rhat, Wgr;
  Array<double,1> Weights, g;
  Vec3 Box, BoxInv;
  void PutInBox(Vec3 &r);
  /// The unnormalized weight function
  inline double phi(double rinv)      { return rinv*rinv*rinv*rinv;}
  /// Its logarithmic derivative
  inline double d_ln_phi(double rinv) { return -4.0*rinv; };
public:
  /// Set the ions and their change
  void Set (const Array<Vec3,1> &rions, const Array<Vec3,1> &delta, Vec3 box);
  Vec3 ForwardWarp(Vec3 r, TinyMatrix<double,3,3> &jmat);
  Vec3 ForwardWarp(Vec3 r);
  /// For a given ending position, this function returns where I must
  /// have warped from.  Solves iteratively.
  Vec3 ReverseWarp(Vec3 rp, TinyMatrix<double,3,3> &jmat);
  Vec3 ReverseWarp(Vec3 rp);
  TinyMatrix<double,3,3> Jmat (Vec3 r);
  /// Returns the inverse of the above
  TinyMatrix<double,3,3> JmatInv (Vec3 r);
};


#endif
