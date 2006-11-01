#ifndef TRICUBIC_B_SPLINE_H
#define TRICUBIC_B_SPLINE_H

#include <blitz/array.h>
#include <blitz/tinymat.h>

using namespace blitz;

class TricubicBspline
{
private:
  TinyMatrix<double,4,4> A, dA, d2A, d3A;
  // The grid sizes
  int Nx, Ny, Nz;
  // The starting and ending values for the uniform grids
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  // The box dimensions and their inverses
  double Lx, LxInv, Ly, LyInv, Lz, LzInv;
  // The grid spacing and inverse
  double dx, dxInv, dy, dyInv, dz, dzInv;
  // The control points
  Array<double,3> P;

  bool Interpolating, Periodic;

  void SolvePeriodicInterp (Array<double,1> data, Array<double,1> p);
  void SolvePeriodicInterp (Array<double,3> &data);
  void MakePeriodic();

  mutable int ix, iy, iz;
  mutable double tx, ty, tz;
  mutable double px[4], py[4], pz[4];
  void Find (double x, double y, double z) const;

public:
  void Init (double xi, double xf, double yi, double yf, double zi, double zf,
	    Array<double,3> &data, bool interp=true, bool periodic=true);
  inline double operator()(double x, double y, double z);
  inline TinyVector<double,3> Grad(double x, double y, double z);
  TricubicBspline();
};

inline void
TricubicBspline::Find(double x, double y, double z) const 
{
  double xDelta = x - xStart;
  double yDelta = y - yStart;
  double zDelta = z - zStart;
  
  if (Periodic) {
//     xDelta -= nearbyint(xDelta*LxInv)*Lx;
//     yDelta -= nearbyint(yDelta*LyInv)*Ly;
//     zDelta -= nearbyint(zDelta*LzInv)*Lz;
    xDelta -= floor(xDelta*LxInv)*Lx;
    yDelta -= floor(yDelta*LyInv)*Ly;
    zDelta -= floor(zDelta*LzInv)*Lz;
  }
  
  double xInt, yInt, zInt;
  tx = modf (xDelta*dxInv, &xInt);
  ty = modf (yDelta*dyInv, &yInt);
  tz = modf (zDelta*dzInv, &zInt);
  ix = (int)xInt;
  iy = (int)yInt;
  iz = (int)zInt;

  px[0] = tx*tx*tx;  py[0] = ty*ty*ty; pz[0] = tz*tz*tz;
  px[1] = tx*tx;     py[1] = ty*ty;    pz[1] = tz*tz;
  px[2] = tx;        py[2] = ty;       pz[2] = tz;
  px[3] = 1.0;       py[3] = 1.0;      pz[3] = 1.0;
}


inline double
TricubicBspline::operator()(double x, double y, double z)
{
  Find(x, y, z);
  double a[4], b[4], c[4];
  a[0] = A(0,0)*px[0]; a[1] =A(0,1)*px[0]; a[2] =A(0,2)*px[0]; a[3] =A(0,3)*px[0];
  a[0]+= A(1,0)*px[1]; a[1]+=A(1,1)*px[1]; a[2]+=A(1,2)*px[1]; a[3]+=A(1,3)*px[1];
  a[0]+= A(2,0)*px[2]; a[1]+=A(2,1)*px[2]; a[2]+=A(2,2)*px[2]; a[3]+=A(2,3)*px[2];
  a[0]+= A(3,0)*px[3]; a[1]+=A(3,1)*px[3]; a[2]+=A(3,2)*px[3]; a[3]+=A(3,3)*px[3];  

  b[0] = A(0,0)*py[0]; b[1] =A(0,1)*py[0]; b[2] =A(0,2)*py[0]; b[3] =A(0,3)*py[0];
  b[0]+= A(1,0)*py[1]; b[1]+=A(1,1)*py[1]; b[2]+=A(1,2)*py[1]; b[3]+=A(1,3)*py[1];
  b[0]+= A(2,0)*py[2]; b[1]+=A(2,1)*py[2]; b[2]+=A(2,2)*py[2]; b[3]+=A(2,3)*py[2];
  b[0]+= A(3,0)*py[3]; b[1]+=A(3,1)*py[3]; b[2]+=A(3,2)*py[3]; b[3]+=A(3,3)*py[3];  

  c[0] = A(0,0)*pz[0]; c[1] =A(0,1)*pz[0]; c[2] =A(0,2)*pz[0]; c[3] =A(0,3)*pz[0];
  c[0]+= A(1,0)*pz[1]; c[1]+=A(1,1)*pz[1]; c[2]+=A(1,2)*pz[1]; c[3]+=A(1,3)*pz[1];
  c[0]+= A(2,0)*pz[2]; c[1]+=A(2,1)*pz[2]; c[2]+=A(2,2)*pz[2]; c[3]+=A(2,3)*pz[2];
  c[0]+= A(3,0)*pz[3]; c[1]+=A(3,1)*pz[3]; c[2]+=A(3,2)*pz[3]; c[3]+=A(3,3)*pz[3];  

  int ix0, ix1, ix2, ix3, iy0, iy1, iy2, iy3, iz0, iz1, iz2, iz3;
  ix0 = ix; ix1=ix+1; ix2=ix+2; ix3=ix+3;
  iy0 = iy; iy1=iy+1; iy2=iy+2; iy3=iy+3;
  iz0 = iz; iz1=iz+1; iz2=iz+2; iz3=iz+3;

  return 
    (a[0]*(b[0]*(c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3))+
	   b[1]*(c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3))+
	   b[2]*(c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3))+
	   b[3]*(c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3)))+
     a[1]*(b[0]*(c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3))+
	   b[1]*(c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3))+
	   b[2]*(c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3))+
	   b[3]*(c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3)))+
     a[2]*(b[0]*(c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3))+
	   b[1]*(c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3))+
	   b[2]*(c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3))+
	   b[3]*(c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3)))+
     a[3]*(b[0]*(c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3))+
	   b[1]*(c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3))+
	   b[2]*(c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3))+
	   b[3]*(c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3))));
}

inline TinyVector<double,3>
TricubicBspline::Grad(double x, double y, double z)
{
  double a[4], b[4], c[4], da[4], db[4], dc[4];
  Find(x, y, z);
  a[0] = A(0,0)*px[0]; a[1] =A(0,1)*px[0]; a[2] =A(0,2)*px[0]; a[3] =A(0,3)*px[0];
  a[0]+= A(1,0)*px[1]; a[1]+=A(1,1)*px[1]; a[2]+=A(1,2)*px[1]; a[3]+=A(1,3)*px[1];
  a[0]+= A(2,0)*px[2]; a[1]+=A(2,1)*px[2]; a[2]+=A(2,2)*px[2]; a[3]+=A(2,3)*px[2];
  a[0]+= A(3,0)*px[3]; a[1]+=A(3,1)*px[3]; a[2]+=A(3,2)*px[3]; a[3]+=A(3,3)*px[3];  

  b[0] = A(0,0)*py[0]; b[1] =A(0,1)*py[0]; b[2] =A(0,2)*py[0]; b[3] =A(0,3)*py[0];
  b[0]+= A(1,0)*py[1]; b[1]+=A(1,1)*py[1]; b[2]+=A(1,2)*py[1]; b[3]+=A(1,3)*py[1];
  b[0]+= A(2,0)*py[2]; b[1]+=A(2,1)*py[2]; b[2]+=A(2,2)*py[2]; b[3]+=A(2,3)*py[2];
  b[0]+= A(3,0)*py[3]; b[1]+=A(3,1)*py[3]; b[2]+=A(3,2)*py[3]; b[3]+=A(3,3)*py[3];  

  c[0] = A(0,0)*pz[0]; c[1] =A(0,1)*pz[0]; c[2] =A(0,2)*pz[0]; c[3] =A(0,3)*pz[0];
  c[0]+= A(1,0)*pz[1]; c[1]+=A(1,1)*pz[1]; c[2]+=A(1,2)*pz[1]; c[3]+=A(1,3)*pz[1];
  c[0]+= A(2,0)*pz[2]; c[1]+=A(2,1)*pz[2]; c[2]+=A(2,2)*pz[2]; c[3]+=A(2,3)*pz[2];
  c[0]+= A(3,0)*pz[3]; c[1]+=A(3,1)*pz[3]; c[2]+=A(3,2)*pz[3]; c[3]+=A(3,3)*pz[3];  

  da[0] = dA(1,0)*px[1]; da[1] =dA(1,1)*px[1]; da[2] =dA(1,2)*px[1]; da[3] =dA(1,3)*px[1];
  da[0]+= dA(2,0)*px[2]; da[1]+=dA(2,1)*px[2]; da[2]+=dA(2,2)*px[2]; da[3]+=dA(2,3)*px[2];
  da[0]+= dA(3,0)*px[3]; da[1]+=dA(3,1)*px[3]; da[2]+=dA(3,2)*px[3]; da[3]+=dA(3,3)*px[3];  

  db[0] = dA(1,0)*py[1]; db[1] =dA(1,1)*py[1]; db[2] =dA(1,2)*py[1]; db[3] =dA(1,3)*py[1];
  db[0]+= dA(2,0)*py[2]; db[1]+=dA(2,1)*py[2]; db[2]+=dA(2,2)*py[2]; db[3]+=dA(2,3)*py[2];
  db[0]+= dA(3,0)*py[3]; db[1]+=dA(3,1)*py[3]; db[2]+=dA(3,2)*py[3]; db[3]+=dA(3,3)*py[3];  

  dc[0] = dA(1,0)*pz[1]; dc[1] =dA(1,1)*pz[1]; dc[2] =dA(1,2)*pz[1]; dc[3] =dA(1,3)*pz[1];
  dc[0]+= dA(2,0)*pz[2]; dc[1]+=dA(2,1)*pz[2]; dc[2]+=dA(2,2)*pz[2]; dc[3]+=dA(2,3)*pz[2];
  dc[0]+= dA(3,0)*pz[3]; dc[1]+=dA(3,1)*pz[3]; dc[2]+=dA(3,2)*pz[3]; dc[3]+=dA(3,3)*pz[3];  

  int ix0, ix1, ix2, ix3, iy0, iy1, iy2, iy3, iz0, iz1, iz2, iz3;
  ix0 = ix; ix1=ix+1; ix2=ix+2; ix3=ix+3;
  iy0 = iy; iy1=iy+1; iy2=iy+2; iy3=iy+3;
  iz0 = iz; iz1=iz+1; iz2=iz+2; iz3=iz+3;


  // Save some operations by factorizing computation.
  TinyMatrix<double,4,4> cP;
  cP(0,0) = c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3);
  cP(0,1) = c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3);
  cP(0,2) = c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3);
  cP(0,3) = c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3);
  cP(1,0) = c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3);
  cP(1,1) = c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3); 
  cP(1,2) = c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3);
  cP(1,3) = c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3);
  cP(2,0) = c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3);
  cP(2,1) = c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3);
  cP(2,2) = c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3);
  cP(2,3) = c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3);
  cP(3,0) = c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3);
  cP(3,1) = c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3);
  cP(3,2) = c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3);
  cP(3,3) = c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3);


  TinyVector<double,3> g;

  g[0] = dxInv *
    (da[0]*(cP(0,0)*b[0]+cP(0,1)*b[1]+cP(0,2)*b[2]+cP(0,3)*b[3]) +
     da[1]*(cP(1,0)*b[0]+cP(1,1)*b[1]+cP(1,2)*b[2]+cP(1,3)*b[3]) +
     da[2]*(cP(2,0)*b[0]+cP(2,1)*b[1]+cP(2,2)*b[2]+cP(2,3)*b[3]) +
     da[3]*(cP(3,0)*b[0]+cP(3,1)*b[1]+cP(3,2)*b[2]+cP(3,3)*b[3]));

  g[1] = dyInv *
    (a[0]*(cP(0,0)*db[0]+cP(0,1)*db[1]+cP(0,2)*db[2]+cP(0,3)*db[3]) +
     a[1]*(cP(1,0)*db[0]+cP(1,1)*db[1]+cP(1,2)*db[2]+cP(1,3)*db[3]) +
     a[2]*(cP(2,0)*db[0]+cP(2,1)*db[1]+cP(2,2)*db[2]+cP(2,3)*db[3]) +
     a[3]*(cP(3,0)*db[0]+cP(3,1)*db[1]+cP(3,2)*db[2]+cP(3,3)*db[3]));

//   g[0] = dxInv *
//     (da[0]*(b[0]*(c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3))+
// 	    b[1]*(c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3))+
// 	    b[2]*(c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3))+
// 	    b[3]*(c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3)))+
//      da[1]*(b[0]*(c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3))+
// 	    b[1]*(c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3))+
// 	    b[2]*(c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3))+
// 	    b[3]*(c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3)))+
//      da[2]*(b[0]*(c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3))+
// 	    b[1]*(c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3))+
// 	    b[2]*(c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3))+
// 	    b[3]*(c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3)))+
//      da[3]*(b[0]*(c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3))+
// 	    b[1]*(c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3))+
// 	    b[2]*(c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3))+
// 	    b[3]*(c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3))));
  
//   g[1] = dyInv * 
//     (a[0]*(db[0]*(c[0]*P(ix0,iy0,iz0)+c[1]*P(ix0,iy0,iz1)+c[2]*P(ix0,iy0,iz2)+c[3]*P(ix0,iy0,iz3))+
// 	   db[1]*(c[0]*P(ix0,iy1,iz0)+c[1]*P(ix0,iy1,iz1)+c[2]*P(ix0,iy1,iz2)+c[3]*P(ix0,iy1,iz3))+
// 	   db[2]*(c[0]*P(ix0,iy2,iz0)+c[1]*P(ix0,iy2,iz1)+c[2]*P(ix0,iy2,iz2)+c[3]*P(ix0,iy2,iz3))+
// 	   db[3]*(c[0]*P(ix0,iy3,iz0)+c[1]*P(ix0,iy3,iz1)+c[2]*P(ix0,iy3,iz2)+c[3]*P(ix0,iy3,iz3)))+
//      a[1]*(db[0]*(c[0]*P(ix1,iy0,iz0)+c[1]*P(ix1,iy0,iz1)+c[2]*P(ix1,iy0,iz2)+c[3]*P(ix1,iy0,iz3))+
// 	   db[1]*(c[0]*P(ix1,iy1,iz0)+c[1]*P(ix1,iy1,iz1)+c[2]*P(ix1,iy1,iz2)+c[3]*P(ix1,iy1,iz3))+
// 	   db[2]*(c[0]*P(ix1,iy2,iz0)+c[1]*P(ix1,iy2,iz1)+c[2]*P(ix1,iy2,iz2)+c[3]*P(ix1,iy2,iz3))+
// 	   db[3]*(c[0]*P(ix1,iy3,iz0)+c[1]*P(ix1,iy3,iz1)+c[2]*P(ix1,iy3,iz2)+c[3]*P(ix1,iy3,iz3)))+
//      a[2]*(db[0]*(c[0]*P(ix2,iy0,iz0)+c[1]*P(ix2,iy0,iz1)+c[2]*P(ix2,iy0,iz2)+c[3]*P(ix2,iy0,iz3))+
// 	   db[1]*(c[0]*P(ix2,iy1,iz0)+c[1]*P(ix2,iy1,iz1)+c[2]*P(ix2,iy1,iz2)+c[3]*P(ix2,iy1,iz3))+
// 	   db[2]*(c[0]*P(ix2,iy2,iz0)+c[1]*P(ix2,iy2,iz1)+c[2]*P(ix2,iy2,iz2)+c[3]*P(ix2,iy2,iz3))+
// 	   db[3]*(c[0]*P(ix2,iy3,iz0)+c[1]*P(ix2,iy3,iz1)+c[2]*P(ix2,iy3,iz2)+c[3]*P(ix2,iy3,iz3)))+
//      a[3]*(db[0]*(c[0]*P(ix3,iy0,iz0)+c[1]*P(ix3,iy0,iz1)+c[2]*P(ix3,iy0,iz2)+c[3]*P(ix3,iy0,iz3))+
// 	   db[1]*(c[0]*P(ix3,iy1,iz0)+c[1]*P(ix3,iy1,iz1)+c[2]*P(ix3,iy1,iz2)+c[3]*P(ix3,iy1,iz3))+
// 	   db[2]*(c[0]*P(ix3,iy2,iz0)+c[1]*P(ix3,iy2,iz1)+c[2]*P(ix3,iy2,iz2)+c[3]*P(ix3,iy2,iz3))+
// 	   db[3]*(c[0]*P(ix3,iy3,iz0)+c[1]*P(ix3,iy3,iz1)+c[2]*P(ix3,iy3,iz2)+c[3]*P(ix3,iy3,iz3))));
  
  g[2] = dzInv * 
    (a[0]*(b[0]*(dc[0]*P(ix0,iy0,iz0)+dc[1]*P(ix0,iy0,iz1)+dc[2]*P(ix0,iy0,iz2)+dc[3]*P(ix0,iy0,iz3))+
	   b[1]*(dc[0]*P(ix0,iy1,iz0)+dc[1]*P(ix0,iy1,iz1)+dc[2]*P(ix0,iy1,iz2)+dc[3]*P(ix0,iy1,iz3))+
	   b[2]*(dc[0]*P(ix0,iy2,iz0)+dc[1]*P(ix0,iy2,iz1)+dc[2]*P(ix0,iy2,iz2)+dc[3]*P(ix0,iy2,iz3))+
	   b[3]*(dc[0]*P(ix0,iy3,iz0)+dc[1]*P(ix0,iy3,iz1)+dc[2]*P(ix0,iy3,iz2)+dc[3]*P(ix0,iy3,iz3)))+
     a[1]*(b[0]*(dc[0]*P(ix1,iy0,iz0)+dc[1]*P(ix1,iy0,iz1)+dc[2]*P(ix1,iy0,iz2)+dc[3]*P(ix1,iy0,iz3))+
	   b[1]*(dc[0]*P(ix1,iy1,iz0)+dc[1]*P(ix1,iy1,iz1)+dc[2]*P(ix1,iy1,iz2)+dc[3]*P(ix1,iy1,iz3))+
	   b[2]*(dc[0]*P(ix1,iy2,iz0)+dc[1]*P(ix1,iy2,iz1)+dc[2]*P(ix1,iy2,iz2)+dc[3]*P(ix1,iy2,iz3))+
	   b[3]*(dc[0]*P(ix1,iy3,iz0)+dc[1]*P(ix1,iy3,iz1)+dc[2]*P(ix1,iy3,iz2)+dc[3]*P(ix1,iy3,iz3)))+
     a[2]*(b[0]*(dc[0]*P(ix2,iy0,iz0)+dc[1]*P(ix2,iy0,iz1)+dc[2]*P(ix2,iy0,iz2)+dc[3]*P(ix2,iy0,iz3))+
	   b[1]*(dc[0]*P(ix2,iy1,iz0)+dc[1]*P(ix2,iy1,iz1)+dc[2]*P(ix2,iy1,iz2)+dc[3]*P(ix2,iy1,iz3))+
	   b[2]*(dc[0]*P(ix2,iy2,iz0)+dc[1]*P(ix2,iy2,iz1)+dc[2]*P(ix2,iy2,iz2)+dc[3]*P(ix2,iy2,iz3))+
	   b[3]*(dc[0]*P(ix2,iy3,iz0)+dc[1]*P(ix2,iy3,iz1)+dc[2]*P(ix2,iy3,iz2)+dc[3]*P(ix2,iy3,iz3)))+
     a[3]*(b[0]*(dc[0]*P(ix3,iy0,iz0)+dc[1]*P(ix3,iy0,iz1)+dc[2]*P(ix3,iy0,iz2)+dc[3]*P(ix3,iy0,iz3))+
	   b[1]*(dc[0]*P(ix3,iy1,iz0)+dc[1]*P(ix3,iy1,iz1)+dc[2]*P(ix3,iy1,iz2)+dc[3]*P(ix3,iy1,iz3))+
	   b[2]*(dc[0]*P(ix3,iy2,iz0)+dc[1]*P(ix3,iy2,iz1)+dc[2]*P(ix3,iy2,iz2)+dc[3]*P(ix3,iy2,iz3))+
	   b[3]*(dc[0]*P(ix3,iy3,iz0)+dc[1]*P(ix3,iy3,iz1)+dc[2]*P(ix3,iy3,iz2)+dc[3]*P(ix3,iy3,iz3))));

  return g;
}
  

  

#endif
