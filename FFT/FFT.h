#ifndef FFT_H
#define FFT_H

#include "../Blitz.h"
#include <fftw3.h>
class FFT3D
{
private:
  complex<double> *rData, *kData;
  fftw_plan r2kPlan, k2rPlan;
  bool Allocated;
  double sqrtNinv;
public:
  Array<complex<double>,3> rBox, kBox;

  void resize (int nx, int ny, int nz);
  void r2k();
  void k2r();

  FFT3D() : Allocated(false)
  {
    // Do nothing for now
  }

};

#endif
