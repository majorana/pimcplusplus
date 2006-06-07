#ifndef LOCAL_POT_FFT_H
#define LOCAL_POT_FFT_H

#include "HamiltonianBase.h"

class LocalPotFFTClass : public VionBase
{
private:
  kSpacePH kPH;
  Array<double,1> Vk;
  zVec Vc;
  FFTBox      &cFFT;
  void SetupkPotentials();
  void SetuprPotentials();
public:
  Array<complex<FFT_FLOAT>,3> Vr;
  void Setup();
  void SetIons (const Array<Vec3, 1>& rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);

  void Apply (const zVec &c, zVec &Hc);
  void Apply (const zVec &c, zVec &Hc,
	      Array<double,3> &VHXC);
  LocalPotFFTClass (Potential &ph, GVecsClass &gvecs, FFTBox &fft) :
    VionBase (gvecs), kPH(ph), cFFT(fft)
  {

  }
};



#endif
