#ifndef MY_TRICUBIC_SPLINE_H
#define MY_TRICUBIC_SPLINE_H

#include "Grid.h"
//#include <blitz/array.h>
//using namespace blitz;

/// Each point of F contains:
/// 0) F(x,y,z)
/// 1) dF/dx
/// 2) dF/dy
/// 3) dF/dz
/// 4) d2F/dxdy
/// 5) d2F/dxdz
/// 6) d2F/dydz
/// 7) d3F/dxdydz

class MyTricubicSpline
{
  inline double p1(double t)
  { return ((t-1.0)*(t-1.0)*(1.0+2.0*t)); }
  inline double p2(double t)
  { return (t*t*(3.0-2.0*t)); }
  inline double q1(double t)
  { return (t*(t-1.0)*(t-1.0)); }
  inline double q2(double t)
  { return (t*t*(t-1.0)); }
  inline double dp1(double t)
  { return (6.0*t*(t-1.0)); }
  inline double dq1(double t)
  { return ((t-1.0)*(3.0*t-1.0)); }
  inline double dp2(double t)
  { return (-dp1(t)); }
  inline double dq2 (double t)
  { return (3.0*t*t - 2.0*t); }
  inline double d2p1(double t)
  { return (12.0*t-6.0); }
  inline double d2q1 (double t)
  { return (6.0*t - 4.0); }
  inline double d2p2 (double t)
  { return (-d2p1(t)); }
  inline double d2q2 (double t)
  { return (6.0*t - 2.0); } 

  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int dest);
  void UpdateY (int source, int dest);
  void UpdateZ (int source, int dest);
  bool UpToDate;
public:
  Array<TinyVector<double,8>,3> F;

  int Nx, Ny, Nz;
  Grid *Xgrid, *Ygrid, *Zgrid;
  TinyVector<Grid*,3> Grids;
  void Update();
  inline double operator()(int ix, int iy, int iz) const
  { return (F(ix,iy,ix)[0]); }
  inline double& operator()(int ix, int iy, int iz) 
  { UpToDate=false; return (F(ix,iy,iz)[0]); }
  inline double operator()(double x, double y, double z);

  MyTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid)
  {
    Xgrid = xgrid; Nx = xgrid->NumPoints;
    Ygrid = ygrid; Ny = ygrid->NumPoints;
    Zgrid = zgrid; Nz = zgrid->NumPoints;
    
    F.resize(Nx,Ny,Nz);
    UpToDate = false;
  }
  
  inline void Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
		    const Array<double,3> &init);

  MyTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid,
		   const Array<double,3> &init)
  {
    Init (xgrid, ygrid, zgrid, init);
  }
  MyTricubicSpline() : UpToDate(false) 
  { /* Do nothing. */ }
};



inline void MyTricubicSpline::Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
				    const Array<double,3> &init)
{
  Xgrid = xgrid; Nx = xgrid->NumPoints;
  Ygrid = ygrid; Ny = ygrid->NumPoints;
  Zgrid = zgrid; Nz = zgrid->NumPoints;

  assert (init.extent(0) == Nx);
  assert (init.extent(1) == Ny);
  assert (init.extent(2) == Nz);

  F.resize(Nx,Ny,Nz);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	F(ix,iy,iz)[0] = init(ix,iy,iz);
  UpToDate = false;
}




inline double MyTricubicSpline::operator() (double x, double y, double z)
{
  if (!UpToDate)
    Update();
  //  Array<double,3> Y(4,4,4);
  TinyMatrix<TinyVector<double,4>,4,4> Y;
  //  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b, c;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);


  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);

  c(0) = p1(w);
  c(1) = p2(w);
  c(2) = l*q1(w);
  c(3) = l*q2(w);
  
  Y(0,0)[0] = F(ix,iy,iz)[0];      //   F
  Y(0,0)[1] = F(ix,iy,iz+1)[0];    //   F
  Y(0,0)[2] = F(ix,iy,iz)[3];      //  dF/dz
  Y(0,0)[3] = F(ix,iy,iz+1)[3];    //  dF/dz
  Y(0,1)[0] = F(ix,iy+1,iz)[0];    //   F
  Y(0,1)[1] = F(ix,iy+1,iz+1)[0];  //   F
  Y(0,1)[2] = F(ix,iy+1,iz)[3];    //  dF/dz
  Y(0,1)[3] = F(ix,iy+1,iz+1)[3];  //  dF/dz
  Y(0,2)[0] = F(ix,iy,iz)[2];      //  dF/dy
  Y(0,2)[1] = F(ix,iy,iz+1)[2];    //  dF/dy
  Y(0,2)[2] = F(ix,iy,iz)[6];      // d2F/dydz
  Y(0,2)[3] = F(ix,iy,iz+1)[6];    // d2F/dydz
  Y(0,3)[0] = F(ix,iy+1,iz)[2];    //  dF/dy
  Y(0,3)[1] = F(ix,iy+1,iz+1)[2];  //  dF/dy
  Y(0,3)[2] = F(ix,iy+1,iz)[6];    // d2F/dydz
  Y(0,3)[3] = F(ix,iy+1,iz+1)[6];  // d2F/dydz

  Y(1,0)[0] = F(ix+1,iy,iz)[0];      //   F
  Y(1,0)[1] = F(ix+1,iy,iz+1)[0];    //   F
  Y(1,0)[2] = F(ix+1,iy,iz)[3];      //  dF/dz
  Y(1,0)[3] = F(ix+1,iy,iz+1)[3];    //  dF/dz
  Y(1,1)[0] = F(ix+1,iy+1,iz)[0];    //   F
  Y(1,1)[1] = F(ix+1,iy+1,iz+1)[0];  //   F
  Y(1,1)[2] = F(ix+1,iy+1,iz)[3];    //  dF/dz
  Y(1,1)[3] = F(ix+1,iy+1,iz+1)[3];  //  dF/dz
  Y(1,2)[0] = F(ix+1,iy,iz)[2];      //  dF/dy
  Y(1,2)[1] = F(ix+1,iy,iz+1)[2];    //  dF/dy
  Y(1,2)[2] = F(ix+1,iy,iz)[6];      // d2F/dydz
  Y(1,2)[3] = F(ix+1,iy,iz+1)[6];    // d2F/dydz
  Y(1,3)[0] = F(ix+1,iy+1,iz)[2];    //  dF/dy
  Y(1,3)[1] = F(ix+1,iy+1,iz+1)[2];  //  dF/dy
  Y(1,3)[2] = F(ix+1,iy+1,iz)[6];    // d2F/dydz
  Y(1,3)[3] = F(ix+1,iy+1,iz+1)[6];  // d2F/dydz

  Y(2,0)[0] = F(ix,iy,iz)[1];      //  dF/dx
  Y(2,0)[1] = F(ix,iy,iz+1)[1];    //  dF/dx
  Y(2,0)[2] = F(ix,iy,iz)[5];      // d2F/dxdz
  Y(2,0)[3] = F(ix,iy,iz+1)[5];    // d2F/dxdz
  Y(2,1)[0] = F(ix,iy+1,iz)[1];    //  dF/dx
  Y(2,1)[1] = F(ix,iy+1,iz+1)[1];  //  dF/dx
  Y(2,1)[2] = F(ix,iy+1,iz)[5];    // d2F/dxdz
  Y(2,1)[3] = F(ix,iy+1,iz+1)[5];  // d2F/dxdz
  Y(2,2)[0] = F(ix,iy,iz)[4];      // d2F/dxdy
  Y(2,2)[1] = F(ix,iy,iz+1)[4];    // d2F/dxdy
  Y(2,2)[2] = F(ix,iy,iz)[7];      // d3F/dxdydz
  Y(2,2)[3] = F(ix,iy,iz+1)[7];    // d3F/dxdydz
  Y(2,3)[0] = F(ix,iy+1,iz)[4];    // d2F/dxdy
  Y(2,3)[1] = F(ix,iy+1,iz+1)[4];  // d2F/dxdy
  Y(2,3)[2] = F(ix,iy+1,iz)[7];    // d3F/dxdydz
  Y(2,3)[3] = F(ix,iy+1,iz+1)[7];  // d3F/dxdydz

  Y(3,0)[0] = F(ix+1,iy,iz)[1];      //  dF/dx
  Y(3,0)[1] = F(ix+1,iy,iz+1)[1];    //  dF/dx
  Y(3,0)[2] = F(ix+1,iy,iz)[5];      // d2F/dxdz
  Y(3,0)[3] = F(ix+1,iy,iz+1)[5];    // d2F/dxdz
  Y(3,1)[0] = F(ix+1,iy+1,iz)[1];    //  dF/dx
  Y(3,1)[1] = F(ix+1,iy+1,iz+1)[1];  //  dF/dx
  Y(3,1)[2] = F(ix+1,iy+1,iz)[5];    // d2F/dxdz
  Y(3,1)[3] = F(ix+1,iy+1,iz+1)[5];  // d2F/dxdz
  Y(3,2)[0] = F(ix+1,iy,iz)[4];      // d2F/dxdy
  Y(3,2)[1] = F(ix+1,iy,iz+1)[4];    // d2F/dxdy
  Y(3,2)[2] = F(ix+1,iy,iz)[7];      // d3F/dxdydz
  Y(3,2)[3] = F(ix+1,iy,iz+1)[7];    // d3F/dxdydz
  Y(3,3)[0] = F(ix+1,iy+1,iz)[4];    // d2F/dxdy
  Y(3,3)[1] = F(ix+1,iy+1,iz+1)[4];  // d2F/dxdy
  Y(3,3)[2] = F(ix+1,iy+1,iz)[7];    // d3F/dxdydz
  Y(3,3)[3] = F(ix+1,iy+1,iz+1)[7];  // d3F/dxdydz
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    Zb_m = 
      b[0]*(Y(m,0)[0]*c[0]+Y(m,0)[1]*c[1]+Y(m,0)[2]*c[2]+Y(m,0)[3]*c[3]) +
      b[1]*(Y(m,1)[0]*c[0]+Y(m,1)[1]*c[1]+Y(m,1)[2]*c[2]+Y(m,1)[3]*c[3]) +
      b[2]*(Y(m,2)[0]*c[0]+Y(m,2)[1]*c[1]+Y(m,2)[2]*c[2]+Y(m,2)[3]*c[3]) +
      b[3]*(Y(m,3)[0]*c[0]+Y(m,3)[1]*c[1]+Y(m,3)[2]*c[2]+Y(m,3)[3]*c[3]);
    val += a(m)*Zb_m;
  }
  

//   double val = 0.0;
//   for (int m=0; m<4; m++) {
//     double Zb_m = 0.0;
//     for (int n=0; n<4; n++)
//       Zb_m += Z(m,n) * b(n);
//     val += Zb_m * a(m);
//   }
  return (val);
}



#endif
