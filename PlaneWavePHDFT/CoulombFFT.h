#ifndef COULOMB_FFT_H
#define COULOMB_FFT_H

#include "HamiltonianBase.h"

class CoulombFFTClass : public VionBase
{
private:
  double Z;
  Array<complex<FFT_FLOAT>,3> Vr;

  void Setup();
  void SetVr();
  FFTBox &FFT;
public:
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);

  CoulombFFTClass (double z, GVecsClass &gvecs, FFTBox &fft) : 
    VionBase (gvecs), Z(z), FFT(fft)
  {
    // nothing for now
  }
};

#endif
