#include "FFT.h"

void TestFFT3D()
{
  FFT3D fft;
  fft.resize(2,2,2);
  fft.kBox(0,0,0) = 1.0;  fft.kBox(0,0,1) = 2.0;
  fft.kBox(0,1,0) = 3.0;  fft.kBox(0,1,1) = 4.0;
  fft.kBox(1,0,0) = 5.0;  fft.kBox(1,0,1) = 6.0;
  fft.kBox(1,1,0) = 7.0;  fft.kBox(1,1,1) = 8.0;

  cerr << "Orignal:\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	fprintf (stderr, "%1.5f\n", fft.kBox(i,j,k).real());
  fft.k2r();
  fft.kBox = 0.0;
  fft.r2k();
    cerr << "FFT^-1(FFT(Original)):\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	fprintf (stderr, "%1.5f\n", fft.kBox(i,j,k).real()/8.0);
}

void TestFFTVec3D()
{
  FFTVec3D fft;
  fft.resize(2,2,2);
  fft.kBox(0,0,0)[0]=1.0;  fft.kBox(0,0,0)[1]=1.1;  fft.kBox(0,0,0)[2]=1.2;  
  fft.kBox(0,0,1)[0]=2.0;  fft.kBox(0,0,1)[1]=2.1;  fft.kBox(0,0,1)[2]=2.2;  
  fft.kBox(0,1,0)[0]=3.0;  fft.kBox(0,1,0)[1]=3.1;  fft.kBox(0,1,0)[2]=3.2;  
  fft.kBox(0,1,1)[0]=4.0;  fft.kBox(0,1,1)[1]=4.1;  fft.kBox(0,1,1)[2]=4.2;  
  fft.kBox(1,0,0)[0]=5.0;  fft.kBox(1,0,0)[1]=5.1;  fft.kBox(1,0,0)[2]=5.2;  
  fft.kBox(1,0,1)[0]=6.0;  fft.kBox(1,0,1)[1]=6.1;  fft.kBox(1,0,1)[2]=6.2;  
  fft.kBox(1,1,0)[0]=7.0;  fft.kBox(1,1,0)[1]=7.1;  fft.kBox(1,1,0)[2]=7.2;  
  fft.kBox(1,1,1)[0]=8.0;  fft.kBox(1,1,1)[1]=8.1;  fft.kBox(1,1,1)[2]=8.2;  

  cerr << "Orignal:\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<3; l++)
	  fprintf (stderr, "%1.5f\n", fft.kBox(i,j,k)[l].real());
  fft.k2r();
  fft.kBox = complex<double>(0.0,0.0);
  fft.r2k();
  cerr << "FFT^-1(FFT(Original)):\n";
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<2; k++)
	for (int l=0; l<3; l++)
	  fprintf (stderr, "%1.5f\n", fft.kBox(i,j,k)[l].real()/8.0);
}


#include <time.h>
#include <stdlib.h>

void TimeTest()
{
  int N = 20;
  FFT3D fft;
  fft.resize(N,N,N);
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
	fft.rBox(i,j,k) = drand48();

  clock_t start,end;
  start = clock();
  for (int n=0; n<1000; n++)
    fft.r2k();
  end = clock();

  fprintf (stderr, "Time for (%d,%d,%d) complex FFT = %1.5f\n",
	   N,N,N, (double)(end-start)/(double)(1000*CLOCKS_PER_SEC));

}


main()
{
  TestFFT3D();
  TestFFTVec3D();
  TimeTest();
}
