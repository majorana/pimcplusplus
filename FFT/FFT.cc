#include "FFT.h"

void 
FFT1D::resize(int n)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*n+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*n+FFT_EXTRA_MEM);
  else
    kData = rData;

  complex<FFT_FLOAT>* rAligned = FFTAlign(rData);
  complex<FFT_FLOAT>* kAligned = FFTAlign(kData);

  Array<complex<FFT_FLOAT>,1> *temp;
  temp = new Array<complex<FFT_FLOAT>,1>(rAligned, shape(n), neverDeleteData);
  rBox.reference (*temp);
  delete temp;
  
  temp = new Array<complex<FFT_FLOAT>,1>(kAligned, shape(n), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = FFTNAME(plan_dft_1d)
    (n, reinterpret_cast<FFTNAME(complex)*>(rAligned), 
     reinterpret_cast<FFTNAME(complex)*>(kAligned), 1, FFTW_MEASURE);
  k2rPlan = FFTNAME(plan_dft_1d) 
    (n, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
     reinterpret_cast<FFTNAME(complex)*>(rAligned), -1,FFTW_MEASURE);

  sqrtNinv = sqrt(1.0/(FFT_FLOAT)n);
  Allocated = true;
}

void 
FFT1D::r2k()
{
  FFTNAME(execute)(r2kPlan);
}


void 
FFT1D::k2r()
{
  FFTNAME(execute)(k2rPlan);
}




void 
FFT3D::resize(int nx, int ny, int nz)
{
  cerr << "FFT box size is " << nx << "x" << ny << "x" << nz << ".\n";
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = (complex<FFT_FLOAT>*) FFT_MALLOC(sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  else
    kData = rData;
  
  complex<FFT_FLOAT>* rAligned = FFTAlign(rData);
  complex<FFT_FLOAT>* kAligned = FFTAlign(kData);
  


  Array<complex<FFT_FLOAT>,3> *temp;
  temp = new Array<complex<FFT_FLOAT>,3>(rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<complex<FFT_FLOAT>,3>(kAligned, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = FFTNAME(plan_dft_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(rAligned), 
     reinterpret_cast<FFTNAME(complex)*>(kAligned), 1, FFTW_MEASURE);
  k2rPlan =  FFTNAME(plan_dft_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
     reinterpret_cast<FFTNAME(complex)*>(rAligned), -1,FFTW_MEASURE);
  
  sqrtNinv = sqrt(1.0/(FFT_FLOAT)(nx*ny*nz));
  Allocated = true;
}


void 
FFT3D::r2k()
{
  FFTNAME(execute)(r2kPlan);
//   for (int i=0; i<kBox.size(); i++)
//     *(kBox.data()+i) *= sqrtNinv;
}


void 
FFT3D::k2r()
{
  FFTNAME(execute)(k2rPlan);
//   for (int i=0; i<rBox.size(); i++)
//     *(rBox.data()+i) *= sqrtNinv;
}


void 
FFTVec3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = ( TinyVector<complex<FFT_FLOAT>,3>*) FFT_MALLOC(3*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = ( TinyVector<complex<FFT_FLOAT>,3>*) FFT_MALLOC(3*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  else
    kData = rData;

  TinyVector<complex<FFT_FLOAT>,3>* rAligned = 
    (TinyVector<complex<FFT_FLOAT>,3>*) FFTAlign((complex<FFT_FLOAT>*) rData);
  TinyVector<complex<FFT_FLOAT>,3>* kAligned = 
    (TinyVector<complex<FFT_FLOAT>,3>*) FFTAlign((complex<FFT_FLOAT>*) kData);

  Array< TinyVector<complex<FFT_FLOAT>,3>,3> *temp;
  temp = new Array< TinyVector<complex<FFT_FLOAT>,3>,3>(rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array< TinyVector<complex<FFT_FLOAT>,3>,3>(kAligned, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = FFTNAME(plan_many_dft)(3, n, 3, 
				   reinterpret_cast<FFTNAME(complex)*>(rAligned),
				   n, 3, 1, 
				   reinterpret_cast<FFTNAME(complex)*>(kAligned), 
				   n, 3, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    FFTNAME(plan_many_dft)(3, n, 3, reinterpret_cast<FFTNAME(complex)*>(kAligned), 
			   n, 3, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(rAligned), n,
			   3, 1, -1, FFTW_MEASURE);
  assert (k2rPlan != NULL);
  sqrtNinv = sqrt(1.0/(FFT_FLOAT)(nx*ny*nz));
  Allocated = true;
}


void 
FFTVec3D::r2k()
{
  FFTNAME(execute)(r2kPlan);
}


void 
FFTVec3D::k2r()
{
  FFTNAME(execute)(k2rPlan);
}




void 
FFTMat3D::resize(int nx, int ny, int nz)
{
  if (Allocated) {
    FFT_FREE(rData);
    if (!InPlace)
      FFT_FREE(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (TinyMatrix<complex<FFT_FLOAT>,3,3>*) 
    FFT_MALLOC(9*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  if (!InPlace)
    kData = (TinyMatrix<complex<FFT_FLOAT>,3,3>*) 
      FFT_MALLOC(9*sizeof(FFTNAME(complex))*nx*ny*nz+FFT_EXTRA_MEM);
  else
    kData = rData;
  
  TinyMatrix<complex<FFT_FLOAT>,3,3>* rAligned = 
    (TinyMatrix<complex<FFT_FLOAT>,3,3>*) FFTAlign((complex<FFT_FLOAT>*)rData);
  TinyMatrix<complex<FFT_FLOAT>,3,3>* kAligned = 
    (TinyMatrix<complex<FFT_FLOAT>,3,3>*) FFTAlign((complex<FFT_FLOAT>*)kData);

  Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3> *temp;
  temp = new Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3>
    (rAligned, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;
  
  temp = new Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3>
    (kAligned, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = 
    FFTNAME(plan_many_dft)(3, n, 9, 
			   reinterpret_cast<FFTNAME(complex)*>(rAligned),
			   n, 9, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(kAligned), 
			   n, 9, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    FFTNAME(plan_many_dft)(3, n, 9, reinterpret_cast<FFTNAME(complex)*>(kAligned),
			   n, 9, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(rAligned), 
			n, 9, 1, -1, FFTW_MEASURE);
  assert (k2rPlan != NULL);

  sqrtNinv = sqrt(1.0/(FFT_FLOAT)(nx*ny*nz));
  Allocated = true;
}


void 
FFTMat3D::r2k()
{
  FFTNAME(execute)(r2kPlan);
}


void 
FFTMat3D::k2r()
{
  FFTNAME(execute)(k2rPlan);
}

