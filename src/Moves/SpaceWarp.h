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


// inline Mat3 operator*(const Mat3 &a, const Mat3 &b)
// {
//   Mat3 c;
//   c(0,0) = a(0,0)*b(0,0)+a(0,1)*b(1,0)+a(0,2)*b(2,0);
//   c(0,1) = a(0,0)*b(0,1)+a(0,1)*b(1,1)+a(0,2)*b(2,1);
//   c(0,2) = a(0,0)*b(0,2)+a(0,1)*b(1,2)+a(0,2)*b(2,2);
//   c(1,0) = a(1,0)*b(0,0)+a(1,1)*b(1,0)+a(1,2)*b(2,0);
//   c(1,1) = a(1,0)*b(0,1)+a(1,1)*b(1,1)+a(1,2)*b(2,1);
//   c(1,2) = a(1,0)*b(0,2)+a(2,1)*b(1,2)+a(1,2)*b(2,2);
//   c(2,0) = a(2,0)*b(0,0)+a(2,1)*b(1,0)+a(2,2)*b(2,0);
//   c(2,1) = a(2,0)*b(0,1)+a(2,1)*b(1,1)+a(2,2)*b(2,1);
//   c(2,2) = a(2,0)*b(0,2)+a(2,1)*b(1,2)+a(2,2)*b(2,2);
//   return c;
// }


// inline double det(const TinyMatrix<double,3,3> &A)
// {
//   return (A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1)) -
// 	  A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0)) +
// 	  A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0)));
// }

// inline Mat3 Inverse (const Mat3 &A)
// {
//   Mat3 Ainv;
//   double dinv = 1.0/det (A);
//   Ainv(0,0) =  dinv*(A(1,1)*A(2,2)-A(1,2)*A(2,1));
//   Ainv(1,0) = -dinv*(A(1,0)*A(2,2)-A(1,2)*A(2,0));
//   Ainv(2,0) =  dinv*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
//   Ainv(0,1) = -dinv*(A(0,1)*A(2,2)-A(0,2)*A(2,1));
//   Ainv(1,1) =  dinv*(A(0,0)*A(2,2)-A(0,2)*A(2,0));
//   Ainv(2,1) = -dinv*(A(0,0)*A(2,1)-A(0,1)*A(2,0));
//   Ainv(0,2) =  dinv*(A(0,1)*A(1,2)-A(0,2)*A(1,1));
//   Ainv(1,2) = -dinv*(A(0,0)*A(1,2)-A(0,2)*A(1,0));
//   Ainv(2,2) =  dinv*(A(0,0)*A(1,1)-A(0,1)*A(1,0));
//   return Ainv;
// }


#endif
