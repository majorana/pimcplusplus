#ifndef PLANE_WAVES_H
#define PLANE_WAVES_H

#include "ConjGrad2.h"

class SystemClass
{
protected:
  GVecsClass GVecs;
  Array<Vec3,1> Rions;
  HamiltonianClass H;
  Array<complex<double>,2> Bands;
  ConjGrad CG;
  Vec3 Box;
  int NumBands;
public:
  void Setup(Vec3 box, double kcut, Potential &ph, bool useFFT=true);
  void Setup(Vec3 box, double kcut, double z, bool useFFT=true);
  void SetIons (const Array<Vec3,1> &rions);
  void DiagonalizeH();

  SystemClass(int numBands) 
    : CG(H, Bands), H(GVecs), NumBands(numBands)
  {

  }
};

#endif
