#ifndef FFT_H
#define FFT_H

#include "../Blitz.h"
#include <fftw3.h>

typedef TinyVector<complex<double>,3>   cVec3;
typedef TinyMatrix<complex<double>,3,3> cMat3;


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

template<int DIM>
class FFTVec1D
{
private:
  TinyVector<complex<double>,DIM> *rData, *kData;
  fftw_plan r2kPlan, k2rPlan;
  bool Allocated, InPlace;
  double sqrtNinv;
public:
  Array<TinyVector<complex<double>,DIM>,1> rBox, kBox;

  void resize (int n)
  {
    if (Allocated) {
      fftw_free(rData);
      if (!InPlace)
	fftw_free(kData);
      fftw_destroy_plan(r2kPlan);
      fftw_destroy_plan(k2rPlan);
    }
    rData = (TinyVector<complex<double>,DIM>*) 
      fftw_malloc(DIM*sizeof(fftw_complex)*n);
    if (!InPlace)
      kData = (TinyVector<complex<double>,DIM>*) 
	fftw_malloc(DIM*sizeof(fftw_complex)*n);
    else
      kData = rData;
    
    Array<TinyVector<complex<double>,DIM>,1> *temp;
    temp = new Array<TinyVector<complex<double>,DIM>,1>
      (rData, shape(n), neverDeleteData);
    rBox.reference (*temp);
    delete temp;
    
    temp = new Array<TinyVector<complex<double>,DIM>,1>
      (kData, shape(n), neverDeleteData);
    kBox.reference(*temp);
    delete temp;
    
    r2kPlan =      
      fftw_plan_many_dft (1, &n, DIM, reinterpret_cast<fftw_complex*>(rData),
			  &n, DIM, 1,  reinterpret_cast<fftw_complex*>(kData), 
			  &n, DIM, 1, 1, FFTW_MEASURE);
    assert (r2kPlan != NULL);
    k2rPlan = 
      fftw_plan_many_dft (1, &n, DIM, reinterpret_cast<fftw_complex*>(kData),
			  &n, DIM, 1, reinterpret_cast<fftw_complex*>(rData), 
			  &n, DIM, 1, -1, FFTW_MEASURE);
    assert (k2rPlan != NULL);
    
    sqrtNinv = sqrt(1.0/(double)n);
    Allocated = true;
  }

  inline int size()
  { return rBox.size(); }
  void r2k() { fftw_execute(r2kPlan); }
  void k2r() { fftw_execute(k2rPlan); }

  FFTVec1D(bool inPlace=true) : Allocated(false), InPlace(inPlace)
  {
    // Do nothing for now
  }
  ~FFTVec1D()
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


class FFT3Df
{
private:
  complex<float> *rData, *kData;
  fftwf_plan r2kPlan, k2rPlan;
  bool Allocated, InPlace;
  float sqrtNinv;
public:
  Array<complex<float>,3> rBox, kBox;

  void resize (int nx, int ny, int nz);
  inline int size()
  { return rBox.size(); }
  void r2k();
  void k2r();

  FFT3Df(bool inPlace=true) : Allocated(false), InPlace(inPlace)
  {
    // Do nothing for now
  }
  ~FFT3Df()
  {
    if (Allocated) {
      fftwf_free(rData);
      if (!InPlace)
	fftw_free(kData);
      fftwf_destroy_plan(r2kPlan);
      fftwf_destroy_plan(k2rPlan);
    }
  }
};


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
