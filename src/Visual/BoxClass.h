#ifndef BOX_CLASS_H
#define BOX_CLASS_H

#include "OnePath.h"

class BoxClass
{
private:
  Vec3 Box, BoxInv;
  inline void PutInBox(Vec3 &r, int dim);
  bool BreakSegment (Vec3 &r1, Vec3 &r2, Vec3 &wall1, Vec3 &wall2);
  bool BreakSegment (Vec3 &r1, Vec3 &r2, Vec3 &wall1, Vec3 &wall2, int dim);
public:
  inline void PutInBox(Vec3 &r);

  inline void Set (Vec3 box) 
  { 
    Box = box; 
    BoxInv[0]=1.0/Box[0]; BoxInv[1]=1.0/Box[1]; BoxInv[2]=1.0/Box[2];
  }
  inline void Set (double lx, double ly, double lz)
  {
    Box[0] = lx; Box[1] = ly; Box[2] = lz;
    BoxInv[0]=1.0/lx; BoxInv[1]=1.0/ly; BoxInv[2]=1.0/lz;
  }
  inline double operator[](int i) const
  { return Box[i]; }
  inline operator Vec3() const
  { return Box; }

  void PutPathsInBox (vector<OnePath*>& inList);
  void PutInBox (Array<double,1> r);
};

inline void BoxClass::PutInBox (Vec3 &r)
{
  for (int i=0; i<3; i++) {
    double n = -floor(r[i]*BoxInv[i]+0.5);
    r[i] += n*Box[i];
  }
}

inline 
void BoxClass::PutInBox (Array<double,1> r)
{
  for (int i=0; i<3; i++) {
    double n = -floor(r(i)*BoxInv[i]+0.5);
    r(i) += n*Box[i];
  }
}

inline void BoxClass::PutInBox (Vec3 &r, int dim)
{
  double n = -floor(r[dim]*BoxInv[dim]+0.5);
  r[dim] += n*Box[dim];
}


#endif
