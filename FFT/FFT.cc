#include "FFT.h"

void FFT3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    fftw_free(rData);
    fftw_free(kData);
    fftw_destroy_plan(r2kPlan);
    fftw_destroy_plan(k2rPlan);
  }
  rData = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
  kData = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);

  Array<complex<double>,3> *temp;
  temp = new Array<complex<double>,3>(rData, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<complex<double>,3>(kData, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = 
    fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(rData), 
		     reinterpret_cast<fftw_complex*>(kData), 1, FFTW_MEASURE);
  k2rPlan = 
    fftw_plan_dft_3d (nx, ny, nz, reinterpret_cast<fftw_complex*>(kData), 
		      reinterpret_cast<fftw_complex*>(rData), -1,FFTW_MEASURE);

  sqrtNinv = sqrt(1.0/(double)(nx*ny*nz));
  Allocated = true;
}


void FFT3D::r2k()
{
  fftw_execute(r2kPlan);
  for (int i=0; i<kBox.size(); i++)
    *(kBox.data()+i) *= sqrtNinv;
}


void FFT3D::k2r()
{
  fftw_execute(k2rPlan);
  for (int i=0; i<rBox.size(); i++)
    *(rBox.data()+i) *= sqrtNinv;
}

