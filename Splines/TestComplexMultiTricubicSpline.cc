/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "MyTricubicSpline.h"
#include "ComplexMultiTricubicSpline.h"

void ValTest()
{
  int N = 10;
  const int numSplines = 4;
  Array<MyTricubicSpline,1> MySplines(2*numSplines);
  ComplexMultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(0.0, 4.0, N);
  yGrid.Init(0.0, 5.0, N);
  zGrid.Init(0.0, 6.0, N);
  
  Array<complex<double>,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N; ix++)
    for (int iy=0; iy<N; iy++)
      for (int iz=0; iz<N; iz++) 
	for (int n=0; n<numSplines; n++)
	  initData(ix,iy,iz,n) = complex<double> (drand48(), drand48());

  cerr << "initData filled.\n";
  for (int n=0; n<numSplines; n++) {
    Array<double,3> temp(N,N,N);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).real();
    MySplines(2*n).Init (&xGrid, &yGrid, &zGrid, temp);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).imag();
    MySplines(2*n+1).Init (&xGrid, &yGrid, &zGrid, temp);
  }
  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData);

  FILE *f1 = fopen ("tri1.dat", "w");
  FILE *f2 = fopen ("tri2.dat", "w");

  Array<complex<double>,1> vals(numSplines);
  for (double u=0.0; u<=1.0; u+=0.001) {
    double x = xGrid.Start + u*(xGrid.End-xGrid.Start);
    double y = xGrid.Start + u*(yGrid.End-yGrid.Start);
    double z = xGrid.Start + u*(zGrid.End-zGrid.Start);
    fprintf (f1, "%1.12f ", u);
    fprintf (f2, "%1.12f ", u);
    MultiSpline(x,y,z,vals);
    for (int n=0; n<numSplines; n++) {
      fprintf (f1, "%1.16e %1.16e ", MySplines(2*n)(x,y,z), MySplines(2*n+1)(x,y,z));
      fprintf (f2, "%1.16e %1.16e ", vals(n).real(), vals(n).imag());
    }
    fprintf (f1, "\n");
    fprintf (f2, "\n");
  }
  fclose(f1);
  fclose(f2);
}


void PeriodicTest()
{
  int N = 10;
  const int numSplines = 4;
  ComplexMultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(0.0, 4.0, N);
  yGrid.Init(0.0, 5.0, N);
  zGrid.Init(0.0, 6.0, N);
  
  Array<complex<double>,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N-1; ix++)
    for (int iy=0; iy<N-1; iy++)
      for (int iz=0; iz<N-1; iz++) 
	for (int n=0; n<numSplines; n++)
	  initData(ix,iy,iz,n) = complex<double>(drand48(), drand48());

  MakePeriodic (initData);

  cerr << "periodic initData filled.\n";

  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData, true);

  FILE *fout = fopen ("MultiPeriodic.dat", "w");

  Array<complex<double>,1> vals(numSplines);
  for (double u=0.0; u<=1.0; u+=0.001) {
    double x = xGrid.Start + u*(xGrid.End-xGrid.Start);
    double y = xGrid.Start + u*(yGrid.End-yGrid.Start);
    double z = xGrid.Start + u*(zGrid.End-zGrid.Start);
    fprintf (fout, "%1.12f ", u);
    MultiSpline(x,y,z,vals);
    for (int n=0; n<numSplines; n++) 
      fprintf (fout, "%1.16e %1.16e", vals(n).real(), vals(n).imag());
    fprintf (fout, "\n");
  }
  fclose(fout);
}




void SpeedTest()
{
  int N = 30;
  const int numSplines = 8;
  Array<MyTricubicSpline,1> MySplines(2*numSplines);
  ComplexMultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(0.0, 4.0, N);
  yGrid.Init(0.0, 5.0, N);
  zGrid.Init(0.0, 6.0, N);
  
  Array<complex<double>,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N; ix++)
    for (int iy=0; iy<N; iy++)
      for (int iz=0; iz<N; iz++) 
	for (int n=0; n<numSplines; n++) 
	  initData(ix,iy,iz,n) = complex<double> (drand48(), drand48());

  for (int n=0; n<numSplines; n++) {
    Array<double,3> temp(N,N,N);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).real();
    MySplines(2*n).Init (&xGrid, &yGrid, &zGrid, temp);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).imag();
    MySplines(2*n+1).Init (&xGrid, &yGrid, &zGrid, temp);
  }

  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData);

  Array<complex<double>,1> vals(numSplines);
  int numEvals = 8*100000;
  clock_t start, end;

  for (int j=0; j<numSplines; j++) 
    vals(j) = complex<double> (MySplines(2*j)(0.1,0.2,0.3), MySplines(2*j+1)(0.1,0.2,0.3));


  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    for (int j=0; j<numSplines; j++)
      vals(j) = complex<double>(MySplines(2*j)(x,y,z), MySplines(2*j+1)(x,y,z));
  }
  end = clock();
  fprintf (stderr, "MySplines time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);

  MultiSpline(0.1, 0.2, 0.3, vals);
  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    MultiSpline(x,y,z,vals);
  }
  end = clock();
  fprintf (stderr, "MultiSpline time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);
}


// extern "C" void z3spline_ (double *x, double *y, double *z,
// 			   double *x0, double *dx, int *nx,
// 			   double *y0, double *dy, int *ny,
// 			   double *z0, double *dz, int *nz,
// 			   void *F, int *num, void *vals);
// extern "C" void z3valgrad_ (double *x, double *y, double *z,
// 			   double *x0, double *dx, int *nx,
// 			   double *y0, double *dy, int *ny,
// 			   double *z0, double *dz, int *nz,
// 			   void *F, int *num, void *vals, void *grads);


void GradValTest()
{
  int N = 100;
  const int numSplines =8;
  Array<MyTricubicSpline,1> MySplines(2*numSplines);
  ComplexMultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(-4.0, 4.0, N);
  yGrid.Init(-5.0, 5.0, N);
  zGrid.Init(-6.0, 6.0, N);
  
  Array<complex<double>,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N; ix++)
    for (int iy=0; iy<N; iy++)
      for (int iz=0; iz<N; iz++) 
	for (int n=0; n<numSplines; n++) 
	  initData(ix,iy,iz,n) = complex<double>(drand48(), drand48());

  for (int n=0; n<numSplines; n++) {
    Array<double,3> temp(N,N,N);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).real();
    MySplines(2*n).Init (&xGrid, &yGrid, &zGrid, temp);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).imag();
    MySplines(2*n+1).Init (&xGrid, &yGrid, &zGrid, temp);
  }

  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData);

  Array<complex<double>,1> vals(numSplines), fvals(numSplines);
  Array<cVec3,1>          grads(numSplines), fgrads(numSplines);
  int numTests = 1000;

  for (int i=0; i<numTests; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    MultiSpline.ValGrad(x,y,z,vals,grads);
    MultiSpline.FValGrad(x,y,z,fvals,fgrads);
    for (int j=0; j<numSplines; j++) {
      complex<double> diff = vals(j) - fvals(j);
      complex<double> v2 = 
	complex<double>(MySplines(2*j)(x,y,z),MySplines(2*j+1)(x,y,z));
      complex<double> diff2 = v2-vals(j);
      if (fabs(real(diff*conj(diff))) > 1.e-14) {
	cerr << "x=" << x << " y=" << y << " z=" << z << endl;
	cerr << "C++ = " << vals(j) << endl;
	cerr << "F77 = " << fvals(j) << endl;
      }
      if (fabs(real(diff*conj(diff2))) > 1.e-14) {
	cerr << "x=" << x << " y=" << y << " z=" << z << endl;
	cerr << "Multi= " << vals(j) << endl;
	cerr << "My   = " << v2 << endl;
      }
    }
    for (int j=0; j<numSplines; j++) 
      for (int k=0; k<3; k++) {
	complex<double> diff = grads(j)[k] - fgrads(j)[k];
	if (fabs(real(diff*conj(diff))) > 1.e-14) {
	  cerr << "x=" << x << " y=" << y << " z=" << z << endl;
	  cerr << "C++ = " << grads(j)[k] << endl;
	  cerr << "F77 = " << fgrads(j)[k] << endl;
	}
	complex<double> v2(MySplines(2*j).Grad(x,y,z)[k],
			   MySplines(2*j+1).Grad(x,y,z)[k]);
	complex<double> diff2 = grads(j)[k] - v2;
	if (fabs(real(diff*conj(diff))) > 1.e-14) {
	  cerr << "x=" << x << " y=" << y << " z=" << z << endl;
	  cerr << "Multi = " << grads(j)[k] << endl;
	  cerr << "My    = " << v2          << endl;
	}
      }
  }
}

void GradSpeedTest()
{
  int N = 30;
  const int numSplines = 8;
  Array<MyTricubicSpline,1> MySplines(2*numSplines);
  ComplexMultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(0.0, 4.0, N);
  yGrid.Init(0.0, 5.0, N);
  zGrid.Init(0.0, 6.0, N);
  
  Array<complex<double>,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N; ix++)
    for (int iy=0; iy<N; iy++)
      for (int iz=0; iz<N; iz++) 
	for (int n=0; n<numSplines; n++) 
	  initData(ix,iy,iz,n) = complex<double>(drand48(), drand48());

  for (int n=0; n<numSplines; n++) {
    Array<double,3> temp(N,N,N);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).real();
    MySplines(2*n).Init (&xGrid, &yGrid, &zGrid, temp);
    for (int ix=0; ix<N; ix++)
      for (int iy=0; iy<N; iy++)
	for (int iz=0; iz<N; iz++)
	  temp(ix,iy,iz) = initData(ix,iy,iz,n).imag();
    MySplines(2*n+1).Init (&xGrid, &yGrid, &zGrid, temp);
  }

  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData);

  Array<complex<double>,1> vals(numSplines), fvals(numSplines);
  Array<cVec3,1>          grads(numSplines), fgrads(numSplines);
  int numEvals = 8*100000;
  clock_t start, end;

  for (int j=0; j<numSplines; j++) 
    vals(j) = complex<double>(MySplines(2*j)(0.1,0.2,0.3), 
			      MySplines(2*j+1)(0.1,0.2,0.3));

  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    for (int j=0; j<numSplines; j++) {
      vals(j)  = 
	complex<double>(MySplines(2*j)(x,y,z),MySplines(2*j+1)(x,y,z));
      Vec3 rGrad, iGrad;
      rGrad = MySplines(2*j).Grad(x,y,z);
      iGrad = MySplines(2*j+1).Grad(x,y,z);
      grads(j)[0] = complex<double>(rGrad[0],iGrad[0]);
      grads(j)[1] = complex<double>(rGrad[0],iGrad[1]);
      grads(j)[2] = complex<double>(rGrad[0],iGrad[2]);
    }
  }
  end = clock();
  fprintf (stderr, "MySplines value and gradient time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);

  MultiSpline(0.1, 0.2, 0.3, vals);


//   for (int i=0; i<N; i++)
//     for (int j=0; j<N; j++)
//       for (int k=0; k<N; k++)
// 	for (int n=0; n<numSplines; n++)
// 	  for (int m=0; m<8; m++) {
// 	    fprintf (stderr, "My   (%d) = (%1.8f, %1.8f)\n", m,
// 		     MySplines(2*n).F(i,j,k)[m],
// 		     MySplines(2*n+1).F(i,j,k)[m]);
// 	    fprintf (stderr, "Multi(%d) = (%1.8f, %1.8f)\n", m,
// 		     MultiSpline.F(i,j,k,n)[m][0],		     
// 		     MultiSpline.F(i,j,k,n)[m][1]);
// 	  }


  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    MultiSpline.ValGrad(x,y,z,vals,grads);
  }
  end = clock();
  fprintf (stderr, "MultiSpline value and gradient time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);

  double x0, y0, z0, dx, dy, dz;
  x0 = xGrid.Start; y0 = yGrid.Start; z0 = zGrid.Start;
  dx = xGrid(1)-xGrid(0); dy=yGrid(1)-yGrid(0); dz=zGrid(1)-zGrid(0);
  int nx,ny,nz;
  nx=N; ny=N; nz=N;
  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
//     z3valgrad_(&x,&y,&z,&x0,&dx,&nx,&y0,&dy,&ny,&z0,&dz,&nz,
// 	       MultiSpline.F.data(), (int *)&numSplines, vals.data(),
// 	       grads.data());
    MultiSpline.FValGrad(x,y,z,vals,grads);
  }
  end = clock();
  fprintf (stderr, "fortran ComplexMultiSpline val + grad time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);

}


main()
{
//   PeriodicTest();
//   ValTest();
   GradValTest();
//  GradSpeedTest();
}
  
