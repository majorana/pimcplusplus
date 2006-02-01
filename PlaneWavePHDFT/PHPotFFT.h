#ifndef PH_POT_FFT_H
#define PH_POT_FFT_H

#include "HamiltonianBase.h"


class PHPotFFTClass : public VionBase
{
private:
  kSpacePH kPH;
  Array<Mat3,1> Fk;
  Array<cMat3,3> Fr;
  Array<double,1> Vk;
  Array<complex<double>,3> Vr;
  zVec Vc;
  zVecVec Gc;
  FFTBox      &cFFT;
  FFTVecBox VecFFT;
  FFTMatBox MatFFT;
  void SetupkPotentials();
  void SetuprPotentials();
public:
  void Setup();
  void SetIons (const Array<Vec3, 1>& rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);

  void Apply (const zVec &c, zVec &Hc);
  PHPotFFTClass (Potential &ph, GVecsClass &gvecs, FFTBox &fft) :
    VionBase (gvecs), kPH(ph), cFFT(fft), VecFFT(gvecs), MatFFT(gvecs)
  {

  }
};

#endif
