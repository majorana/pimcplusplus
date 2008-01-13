#ifndef LATTICE_CLASS_H
#define LATTICE_CLASS_H

#include "../Blitz.h"


class LatticeClass
{
private:
  Mat3 Lattice, LatticeInv;
  // RecipLattice = 2 pi Transpose(Inverse(Lattice))
  Mat3 RecipLattice;
public:
  // Returns i'th lattice vector
  inline Vec3 operator[](int i);
  inline double operator()(int i, int j);
  inline Vec3 a(int i);
  inline Vec3 b(int i);
  inline void Set (Mat3 lattice);
  inline Vec3 Cart2Reduced (Vec3 r);
  inline Vec3 Reduced2Cart (Vec3 r);
  // Returns the radius of the largest inscribed sphere
  inline double rMax();
  inline TinyVector<int,3> MinIndices(double kcut);

  LatticeClass (Mat3 lattice)
  { Set (lattice); }
  LatticeClass()
  { /* do nothing */  }
};

inline Vec3
LatticeClass::a(int i)
{
  return Vec3 (Lattice(i,0), Lattice(i,1), Lattice(i,2));
}

inline Vec3
LatticeClass::b(int i)
{
  return Vec3 (RecipLattice(i,0), RecipLattice(i,1), RecipLattice(i,2));
}

inline Vec3
LatticeClass::operator[](int i)
{
  return Vec3 (Lattice(i,0), Lattice(i,1), Lattice(i,2));
}

inline double
LatticeClass::operator()(int i, int j)
{
  return Lattice(i,j);
}

inline void
LatticeClass::Set(Mat3 lattice)
{
  Lattice = lattice;
  LatticeInv = Inverse(lattice);
  RecipLattice = 2.0*M_PI*Transpose(LatticeInv);
}

inline Vec3
LatticeClass::Cart2Reduced (Vec3 r)
{
  return LatticeInv * r;
}

inline Vec3
LatticeClass::Reduced2Cart (Vec3 u)
{
  return u * Lattice;
}

inline double
LatticeClass::rMax()
{
  double r = 1.0e200;
  Vec3 c = 0.5*((*this)[0] + (*this)[1] + (*this)[2]);
  for (int i=0; i<3; i++) {
    Vec3 v1 = (*this)[(i+1)%3];
    Vec3 v2 = (*this)[(i+2)%3];
    
    double beta1 = (dot(v2,v2)*dot(c,v1) - dot (v1,v2)*dot(c,v2))/
	(dot(v1,v1)*dot(v2,v2) - dot(v1,v2)*dot(v1,v2));
    double beta2 = (dot(v1,v1)*dot(c,v2) - dot (v1,v2)*dot(c,v1))/
	(dot(v1,v1)*dot(v2,v2) - dot(v1,v2)*dot(v1,v2));

    Vec3 p = beta1*v1 + beta2*v2;
    double dist = sqrt(dot(c-p, c-p));
    r = min (r, dist);
  }
  return r;
}


// This returns the triple of integers defining the multiples of the
// reciprocal lattice vectors necessary to contain the G-vectors with
// a given reciprocal-space cutoff
TinyVector<int,3>
LatticeClass::MinIndices(double kcut)
{
  Vec3 I;
  for (int i=0; i<3; i++) {
    Vec3 v1 = b((i+1)%3);
    Vec3 v2 = b((i+2)%3);
    Vec3 perpDir = cross(v1,v2);
    perpDir = (1.0/sqrt(dot(perpDir, perpDir)))*perpDir;
    double perpDist = fabs(dot(b(i), perpDir));
    // We add one just do be safe
    I[i] = (int)ceil(kcut/perpDist)+1;
  }
  return I;
}

#endif
