#include "FFT.h"

void 
FFT1D::resize(int n)
{
  if (Allocated) {
    FFTNAME(free)(rData);
    if (!InPlace)
      FFTNAME(free)(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (complex<FFT_FLOAT>*) FFTNAME(malloc)(sizeof(FFTNAME(complex))*n);
  if (!InPlace)
    kData = (complex<FFT_FLOAT>*) FFTNAME(malloc)(sizeof(FFTNAME(complex))*n);
  else
    kData = rData;

  Array<complex<FFT_FLOAT>,1> *temp;
  temp = new Array<complex<FFT_FLOAT>,1>(rData, shape(n), neverDeleteData);
  rBox.reference (*temp);
  delete temp;
  
  temp = new Array<complex<FFT_FLOAT>,1>(kData, shape(n), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = FFTNAME(plan_dft_1d)
    (n, reinterpret_cast<FFTNAME(complex)*>(rData), 
     reinterpret_cast<FFTNAME(complex)*>(kData), 1, FFTW_MEASURE);
  k2rPlan = FFTNAME(plan_dft_1d) 
    (n, reinterpret_cast<FFTNAME(complex)*>(kData), 
     reinterpret_cast<FFTNAME(complex)*>(rData), -1,FFTW_MEASURE);

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
    FFTNAME(free)(rData);
    if (!InPlace)
      FFTNAME(free)(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (complex<FFT_FLOAT>*) 
    FFTNAME(malloc)(sizeof(FFTNAME(complex))*nx*ny*nz);
  if (!InPlace)
    kData = (complex<FFT_FLOAT>*) FFTNAME(malloc)
      (sizeof(FFTNAME(complex))*nx*ny*nz);
  else
    kData = rData;
  
  Array<complex<FFT_FLOAT>,3> *temp;
  temp = new Array<complex<FFT_FLOAT>,3>(rData, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array<complex<FFT_FLOAT>,3>(kData, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  r2kPlan = FFTNAME(plan_dft_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(rData), 
     reinterpret_cast<FFTNAME(complex)*>(kData), 1, FFTW_MEASURE);
  k2rPlan =  FFTNAME(plan_dft_3d)
    (nx, ny, nz, reinterpret_cast<FFTNAME(complex)*>(kData), 
     reinterpret_cast<FFTNAME(complex)*>(rData), -1,FFTW_MEASURE);
  
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
    FFTNAME(free)(rData);
    if (!InPlace)
      FFTNAME(free)(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = ( TinyVector<complex<FFT_FLOAT>,3>*) FFTNAME(malloc)(3*sizeof(FFTNAME(complex))*nx*ny*nz);
  if (!InPlace)
    kData = ( TinyVector<complex<FFT_FLOAT>,3>*) FFTNAME(malloc)(3*sizeof(FFTNAME(complex))*nx*ny*nz);
  else
    kData = rData;

  Array< TinyVector<complex<FFT_FLOAT>,3>,3> *temp;
  temp = new Array< TinyVector<complex<FFT_FLOAT>,3>,3>(rData, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;

  temp = new Array< TinyVector<complex<FFT_FLOAT>,3>,3>(kData, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = FFTNAME(plan_many_dft)(3, n, 3, 
				   reinterpret_cast<FFTNAME(complex)*>(rData),
				   n, 3, 1, 
				   reinterpret_cast<FFTNAME(complex)*>(kData), 
				   n, 3, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    FFTNAME(plan_many_dft)(3, n, 3, reinterpret_cast<FFTNAME(complex)*>(kData),
			   n, 3, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(rData), n,
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
    FFTNAME(free)(rData);
    if (!InPlace)
      FFTNAME(free)(kData);
    FFTNAME(destroy_plan)(r2kPlan);
    FFTNAME(destroy_plan)(k2rPlan);
  }
  rData = (TinyMatrix<complex<FFT_FLOAT>,3,3>*) 
    FFTNAME(malloc)(9*sizeof(FFTNAME(complex))*nx*ny*nz);
  if (!InPlace)
    kData = (TinyMatrix<complex<FFT_FLOAT>,3,3>*) 
      FFTNAME(malloc)(9*sizeof(FFTNAME(complex))*nx*ny*nz);
  else
    kData = rData;
  
  Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3> *temp;
  temp = new Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3>
    (rData, shape(nx,ny,nz), neverDeleteData);
  rBox.reference (*temp);
  delete temp;
  
  temp = new Array<TinyMatrix<complex<FFT_FLOAT>,3,3>,3>
    (kData, shape(nx,ny,nz), neverDeleteData);
  kBox.reference(*temp);
  delete temp;

  int n[3] = {nx, ny, nz};
  r2kPlan = 
    FFTNAME(plan_many_dft)(3, n, 9, 
			   reinterpret_cast<FFTNAME(complex)*>(rData),
			   n, 9, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(kData), 
			   n, 9, 1, 1, FFTW_MEASURE);
  assert (r2kPlan != NULL);
  k2rPlan = 
    FFTNAME(plan_many_dft)(3, n, 9, reinterpret_cast<FFTNAME(complex)*>(kData),
			   n, 9, 1, 
			   reinterpret_cast<FFTNAME(complex)*>(rData), 
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

