#include "FFT.h"

void 
FFT1D::resize(int n)
{
  if (Allocated) {
    fftw_free(rData);
    if (!InPlace)
      fftw_free(kData);
    fftw_destroy_plan(r2kPlan);
    fftw_destroy_plan(k2rPlan);
  }
  rData = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*n);
  if (!InPlace)
    kData = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*n);
  else
    kData = rData;

  Array<complex<double>,1> *temp;
  temp = new Array<complex<double>,1>(rData, shape(n), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<complex<double>,1>(kData, shape(n), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = 
    fftw_plan_dft_1d(n, reinterpret_cast<fftw_complex*>(rData), 
		     reinterpret_cast<fftw_complex*>(kData), 1, FFTW_MEASURE);
  k2rPlan = 
    fftw_plan_dft_1d (n, reinterpret_cast<fftw_complex*>(kData), 
		      reinterpret_cast<fftw_complex*>(rData), -1,FFTW_MEASURE);

  sqrtNinv = sqrt(1.0/(double)n);
  Allocated = true;
}

void 
FFT1D::r2k()
{
  fftw_execute(r2kPlan);
}


void 
FFT1D::k2r()
{
  fftw_execute(k2rPlan);
}




void 
FFT3D::resize(int nx, int ny, int nz)
{
  cerr << "FFT box size is " << nx << "x" << ny << "x" << nz << ".\n";
  if (Allocated) {
    fftw_free(rData);
    if (!InPlace)
      fftw_free(kData);
    fftw_destroy_plan(r2kPlan);
    fftw_destroy_plan(k2rPlan);
  }
  rData = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
  if (!InPlace)
    kData = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
  else
    kData = rData;

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


void 
FFT3D::r2k()
{
  fftw_execute(r2kPlan);
//   for (int i=0; i<kBox.size(); i++)
//     *(kBox.data()+i) *= sqrtNinv;
}


void 
FFT3D::k2r()
{
  fftw_execute(k2rPlan);
//   for (int i=0; i<rBox.size(); i++)
//     *(rBox.data()+i) *= sqrtNinv;
}

void 
FFTVec3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    fftw_free(rData);
    if (!InPlace)
      fftw_free(kData);
    fftw_destroy_plan(r2kPlan);
    fftw_destroy_plan(k2rPlan);
  }
  rData = (cVec3*) fftw_malloc(3*sizeof(fftw_complex)*nx*ny*nz);
  if (!InPlace)
    kData = (cVec3*) fftw_malloc(3*sizeof(fftw_complex)*nx*ny*nz);
  else
    kData = rData;

  Array<cVec3,3> *temp;
  temp = new Array<cVec3,3>(rData, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<cVec3,3>(kData, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = 
    fftw_plan_many_dft (3, n, 3, 
			reinterpret_cast<fftw_complex*>(rData),
			n, 3, 1, 
			reinterpret_cast<fftw_complex*>(kData), 
			n, 3, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    fftw_plan_many_dft (3, n, 3, reinterpret_cast<fftw_complex*>(kData),
			n, 3, 1, 
			reinterpret_cast<fftw_complex*>(rData), n,
			3, 1, -1, FFTW_MEASURE);
  assert (k2rPlan != NULL);

//   r2kPlan = 
//     fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(rData), 
// 		     reinterpret_cast<fftw_complex*>(kData), 1, FFTW_MEASURE);
//   k2rPlan = 
//     fftw_plan_dft_3d (nx, ny, nz, reinterpret_cast<fftw_complex*>(kData), 
// 		      reinterpret_cast<fftw_complex*>(rData), -1,FFTW_MEASURE);

  sqrtNinv = sqrt(1.0/(double)(nx*ny*nz));
  Allocated = true;
}


void 
FFTVec3D::r2k()
{
  fftw_execute(r2kPlan);
}


void 
FFTVec3D::k2r()
{
  fftw_execute(k2rPlan);
}




void 
FFTMat3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    fftw_free(rData);
    if (!InPlace)
      fftw_free(kData);
    fftw_destroy_plan(r2kPlan);
    fftw_destroy_plan(k2rPlan);
  }
  rData = (cMat3*) fftw_malloc(9*sizeof(fftw_complex)*nx*ny*nz);
  if (!InPlace)
    kData = (cMat3*) fftw_malloc(9*sizeof(fftw_complex)*nx*ny*nz);
  else
    kData = rData;

  Array<cMat3,3> *temp;
  temp = new Array<cMat3,3>(rData, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<cMat3,3>(kData, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = 
    fftw_plan_many_dft (3, n, 9, 
			reinterpret_cast<fftw_complex*>(rData),
			n, 9, 1, 
			reinterpret_cast<fftw_complex*>(kData), 
			n, 9, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    fftw_plan_many_dft (3, n, 9, reinterpret_cast<fftw_complex*>(kData),
			n, 9, 1, 
			reinterpret_cast<fftw_complex*>(rData), 
			n, 9, 1, -1, FFTW_MEASURE);
  assert (k2rPlan != NULL);

  sqrtNinv = sqrt(1.0/(double)(nx*ny*nz));
  Allocated = true;
}


void 
FFTMat3D::r2k()
{
  fftw_execute(r2kPlan);
}


void 
FFTMat3D::k2r()
{
  fftw_execute(k2rPlan);
}

