#include "TricubicBspline.h"
#include "../PlaneWavePHDFT/FFTBox.h"

void MakeBlipCoefs (FFTBox &fft,  const zVec &c,
		    Array<complex<double>,3> &coefs)
{
  zVec gc(c.size());
  int nx, ny, nz;
  fft.GetDims(nx,ny,nz);
  Vec3 box = fft.GVecs.GetBox();
  Vec3 a;
  a[0] = box[0]/nx;
  a[1] = box[1]/ny;
  a[2] = box[2]/nz;
  
  for (int i=0; i<c.size(); i++) {
    Vec3 G = fft.GVecs(i);
    G[0] *= 1.0*a[0];
    G[1] *= 1.0*a[1];
    G[2] *= 1.0*a[2];
    double gamma = 1.0;
    if (fabs(G[0]) > 1.0e-14)
      gamma *= (3.0/(G[0]*G[0]*G[0]*G[0])*(3.0 - 4.0*cos(G[0]) + cos(2.0*G[0])));
    else
      gamma *= 1.5;
    if (fabs(G[1]) > 1.0e-14)
      gamma *= (3.0/(G[1]*G[1]*G[1]*G[1])*(3.0 - 4.0*cos(G[1]) + cos(2.0*G[1])));
    else
      gamma *= 1.5;
    if (fabs(G[2]) > 1.0e-14)
      gamma *= (3.0/(G[2]*G[2]*G[2]*G[2])*(3.0 - 4.0*cos(G[2]) + cos(2.0*G[2])));
    else
      gamma *= 1.5;

    gc(i) = c(i)/gamma;
  }
  fft.PutkVec (gc);
  fft.k2r();
  coefs.resize(nx,ny,nz);
  coefs = fft.rBox;
}

void
MakeInterpCoefs (FFTBox &fft, const zVec &c,
		 Array<complex<double>,3> &coefs)
{
  fft.kBox = complex<double>();
  fft.PutkVec (c);
  fft.k2r();
  Array<complex<double>,3> data;
  int nx, ny, nz;
  fft.GetDims(nx,ny,nz);
  data.resize(nx,ny,nz);
  data = fft.rBox;
  TricubicBspline<complex<double> > spline;
  Vec3 box = fft.GVecs.GetBox();
  spline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
	       data, true, true);
  coefs.resize(nx,ny,nz);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	coefs(ix,iy,iz) = spline.GetControlPoint (ix,iy,iz);

}

complex<double> eval(double x, double y, double z, 
		     zVec &coefs, GVecsClass &gvecs)
{
  complex<double> val = complex<double>();
  for (int i=0; i<coefs.size(); i++) {
    Vec3 G = gvecs(i);
    double phase = -(x*G[0]+ y *G[1] + z*G[2]);
    val += coefs(i)*complex<double> (cos(phase), sin(phase));
  }
  return val;
}


void
TestBlips()
{
  const double kcut = 3.0*M_PI;
  Vec3 box (10.0, 10.0, 10.0);
  GVecsClass GVecs;
  Vec3 k (0.0, 0.0, 0.0);
  
  GVecs.Set(box, k, kcut, 4.0);
  FFTBox fft(GVecs);
  fft.Setup();
  zVec c(GVecs.size());
  for (int i=0; i<GVecs.size(); i++)
    c(i) = complex<double> (-1.0+2.0*drand48(),-1.0+2.0*drand48());

  Array<complex<double>,3> blipCoefs, interpCoefs;
  int nx, ny, nz;
  fft.GetDims(nx, ny, nz);
  blipCoefs.resize(nx, ny, nz);
  interpCoefs.resize(nx, ny, nz);
  MakeBlipCoefs(fft, c, blipCoefs);
   MakeInterpCoefs(fft, c, interpCoefs);
//   fft.PutkVec(c);
//   fft.k2r();
//   interpCoefs = fft.rBox;
  //blipCoefs   = fft.rBox;
//   for (int ix=0; ix<interpCoefs.extent(0); ix++)
//     for (int iy=0; iy<interpCoefs.extent(1); iy++)
//       for (int iz=0; iz<interpCoefs.extent(2); iz++) {
// 	double x = box[0]*(double)ix/interpCoefs.extent(0);
// 	double y = box[1]*(double)iy/interpCoefs.extent(1);
// 	double z = box[2]*(double)iz/interpCoefs.extent(2);
// 	interpCoefs(ix,iy,iz) = eval(x,y,z,c,GVecs);
//       }
  TricubicBspline<complex<double> > blipSpline, interpSpline;
  blipSpline.Init (0.0, box[0], 0.0, box[1], 0.0, box[2],
		   blipCoefs, false, true);
  interpSpline.Init (0.0, box[1], 0.0, box[1], 0.0, box[2],
		     interpCoefs, false, true);
//   for (int ix=0; ix<blipCoefs.extent(0); ix++)
//     for (int iy=0; iy<blipCoefs.extent(1); iy++)
//       for (int iz=0; iz<blipCoefs.extent(2); iz++)
// 	fprintf (stdout, "%22.16e %22.16e\n", blipCoefs(ix,iy,iz).real(),
// 		 interpCoefs(ix,iy,iz).real());
  double y= 2.3;
  double z = 3.8;
  for (double x=0.0; x<=10.0; x+=0.001)
    fprintf (stdout, "%22.16e %22.16e %22.16e %22.16e\n", x, interpSpline(x,y,z).real(),
	     blipSpline(x,y,z).real(), eval(x,y,z,c,GVecs).real());
	     
}


main()
{
  TestBlips();
}
