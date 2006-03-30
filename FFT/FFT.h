#ifndef FFT_H
#define FFT_H

#include "../Blitz.h"
#include <fftw3.h>

class FFT1D
{
private:
  complex<double> *rData, *kData;
  fftw_plan r2kPlan, k2rPlan;
  bool Allocated, InPlace;
  double sqrtNinv;
public:
  Array<complex<double>,1> rBox, kBox;

  void resize (int n);
  inline int size()
  { return rBox.size(); }
  void r2k();
  void k2r();

  FFT1D(bool inPlace=true) : Allocated(false), InPlace(inPlace)
  {
    // Do nothing for now
  }
  ~FFT1D()
  {
    if (Allocated) {
      fftw_free(rData);
      if (!InPlace)
	fftw_free(kData);
      fftw_destroy_plan(r2kPlan);
      fftw_destroy_plan(k2rPlan);
    }
  }
};


class FFT3D
{
private:
  complex<double> *rData, *kData;
  fftw_plan r2kPlan, k2rPlan;
  bool Allocated, InPlace;
  double sqrtNinv;
public:
  Array<complex<double>,3> rBox, kBox;

  void resize (int nx, int ny, int nz);
  inline int size()
  { return rBox.size(); }
  void r2k();
  void k2r();

  FFT3D(bool inPlace=true) : Allocated(false), InPlace(inPlace)
  {
    // Do nothing for now
  }
  ~FFT3D()
  {
    if (Allocated) {
      fftw_free(rData);
      if (!InPlace)
	fftw_free(kData);
      fftw_destroy_plan(r2kPlan);
      fftw_destroy_plan(k2rPlan);
    }
  }
};

typedef TinyVector<complex<double>,3>   cVec3;
typedef TinyMatrix<complex<double>,3,3> cMat3;

class FFTVec3D
{
private:
  cVec3 *rData, *kData;
  fftw_plan r2kPlan, k2rPlan;
  bool Allocated, InPlace;
  double sqrtNinv;
public:
  Array<cVec3,3> rBox, kBox;

  void resize (int nx, int ny, int nz);
  inline int size()
  { return rBox.size(); }
  void r2k();
  void k2r();

  FFTVec3D(bool inPlace=true) : Allocated(false), InPlace(inPlace)
  {
    // Do nothing for now
  }
  ~FFTVec3D()
  {
    if (Allocated) {
      fftw_free(rData);
      if (!InPlace)
	fftw_free(kData);
      fftw_destroy_plan(r2kPlan);
      fftw_destroy_plan(k2rPlan);
    }
  }
};

class FFTMat3D
{
private:
  cMat3 *rData, *kData;
  fftw_plan r2kPlan, k2rPlan;
  bool Allocated, InPlace;
  double sqrtNinv;
public:
  Array<cMat3,3> rBox, kBox;

  void resize (int nx, int ny, int nz);
  inline int size()
  { return rBox.size(); }
  void r2k();
  void k2r();

  FFTMat3D(bool inPlace=true) : Allocated(false), InPlace(inPlace)
  {
    // Do nothing for now
  }
  ~FFTMat3D()
  {
    if (Allocated) {
      cerr << "Deallocating FFTMat3D.\n";
      fftw_free(rData);
      if (!InPlace)
	fftw_free(kData);
      fftw_destroy_plan(r2kPlan);
      fftw_destroy_plan(k2rPlan);
    }
  }
};


#endif
