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
  
  /// Copy constructor
  inline MyTricubicSpline (const MyTricubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline MyTricubicSpline & operator= (MyTricubicSpline &a);
  inline MyTricubicSpline & operator= (MyTricubicSpline a);


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


inline MyTricubicSpline::MyTricubicSpline(const MyTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
}


inline MyTricubicSpline& MyTricubicSpline::operator=(MyTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}

inline MyTricubicSpline& MyTricubicSpline::operator=(MyTricubicSpline a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}



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
  //TinyMatrix<TinyVector<double,4>,4,4> Y;
  //  TinyMatrix<double,4,4> Z;
  //TinyVector<double,4> a, b, c;

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
  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  double& Y000 = F(ix,iy,iz)[0];      //   F
  double& Y001 = F(ix,iy,iz+1)[0];    //   F
  double& Y002 = F(ix,iy,iz)[3];      //  dF/dz
  double& Y003 = F(ix,iy,iz+1)[3];    //  dF/dz
  double& Y010 = F(ix,iy+1,iz)[0];    //   F
  double& Y011 = F(ix,iy+1,iz+1)[0];  //   F
  double& Y012 = F(ix,iy+1,iz)[3];    //  dF/dz
  double& Y013 = F(ix,iy+1,iz+1)[3];  //  dF/dz
  double& Y020 = F(ix,iy,iz)[2];      //  dF/dy
  double& Y021 = F(ix,iy,iz+1)[2];    //  dF/dy
  double& Y022 = F(ix,iy,iz)[6];      // d2F/dydz
  double& Y023 = F(ix,iy,iz+1)[6];    // d2F/dydz
  double& Y030 = F(ix,iy+1,iz)[2];    //  dF/dy
  double& Y031 = F(ix,iy+1,iz+1)[2];  //  dF/dy
  double& Y032 = F(ix,iy+1,iz)[6];    // d2F/dydz
  double& Y033 = F(ix,iy+1,iz+1)[6];  // d2F/dydz

  double& Y100 = F(ix+1,iy,iz)[0];      //   F
  double& Y101 = F(ix+1,iy,iz+1)[0];    //   F
  double& Y102 = F(ix+1,iy,iz)[3];      //  dF/dz
  double& Y103 = F(ix+1,iy,iz+1)[3];    //  dF/dz
  double& Y110 = F(ix+1,iy+1,iz)[0];    //   F
  double& Y111 = F(ix+1,iy+1,iz+1)[0];  //   F
  double& Y112 = F(ix+1,iy+1,iz)[3];    //  dF/dz
  double& Y113 = F(ix+1,iy+1,iz+1)[3];  //  dF/dz
  double& Y120 = F(ix+1,iy,iz)[2];      //  dF/dy
  double& Y121 = F(ix+1,iy,iz+1)[2];    //  dF/dy
  double& Y122 = F(ix+1,iy,iz)[6];      // d2F/dydz
  double& Y123 = F(ix+1,iy,iz+1)[6];    // d2F/dydz
  double& Y130 = F(ix+1,iy+1,iz)[2];    //  dF/dy
  double& Y131 = F(ix+1,iy+1,iz+1)[2];  //  dF/dy
  double& Y132 = F(ix+1,iy+1,iz)[6];    // d2F/dydz
  double& Y133 = F(ix+1,iy+1,iz+1)[6];  // d2F/dydz

  double& Y200 = F(ix,iy,iz)[1];      //  dF/dx
  double& Y201 = F(ix,iy,iz+1)[1];    //  dF/dx
  double& Y202 = F(ix,iy,iz)[5];      // d2F/dxdz
  double& Y203 = F(ix,iy,iz+1)[5];    // d2F/dxdz
  double& Y210 = F(ix,iy+1,iz)[1];    //  dF/dx
  double& Y211 = F(ix,iy+1,iz+1)[1];  //  dF/dx
  double& Y212 = F(ix,iy+1,iz)[5];    // d2F/dxdz
  double& Y213 = F(ix,iy+1,iz+1)[5];  // d2F/dxdz
  double& Y220 = F(ix,iy,iz)[4];      // d2F/dxdy
  double& Y221 = F(ix,iy,iz+1)[4];    // d2F/dxdy
  double& Y222 = F(ix,iy,iz)[7];      // d3F/dxdydz
  double& Y223 = F(ix,iy,iz+1)[7];    // d3F/dxdydz
  double& Y230 = F(ix,iy+1,iz)[4];    // d2F/dxdy
  double& Y231 = F(ix,iy+1,iz+1)[4];  // d2F/dxdy
  double& Y232 = F(ix,iy+1,iz)[7];    // d3F/dxdydz
  double& Y233 = F(ix,iy+1,iz+1)[7];  // d3F/dxdydz

  double& Y300 = F(ix+1,iy,iz)[1];      //  dF/dx
  double& Y301 = F(ix+1,iy,iz+1)[1];    //  dF/dx
  double& Y302 = F(ix+1,iy,iz)[5];      // d2F/dxdz
  double& Y303 = F(ix+1,iy,iz+1)[5];    // d2F/dxdz
  double& Y310 = F(ix+1,iy+1,iz)[1];    //  dF/dx
  double& Y311 = F(ix+1,iy+1,iz+1)[1];  //  dF/dx
  double& Y312 = F(ix+1,iy+1,iz)[5];    // d2F/dxdz
  double& Y313 = F(ix+1,iy+1,iz+1)[5];  // d2F/dxdz
  double& Y320 = F(ix+1,iy,iz)[4];      // d2F/dxdy
  double& Y321 = F(ix+1,iy,iz+1)[4];    // d2F/dxdy
  double& Y322 = F(ix+1,iy,iz)[7];      // d3F/dxdydz
  double& Y323 = F(ix+1,iy,iz+1)[7];    // d3F/dxdydz
  double& Y330 = F(ix+1,iy+1,iz)[4];    // d2F/dxdy
  double& Y331 = F(ix+1,iy+1,iz+1)[4];  // d2F/dxdy
  double& Y332 = F(ix+1,iy+1,iz)[7];    // d3F/dxdydz
  double& Y333 = F(ix+1,iy+1,iz+1)[7];  // d3F/dxdydz
  
  double val = 
    a0*
    (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
     b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
     b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
     b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3)) +
    a1 *
    (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
     b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
     b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
     b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
    a2 *
    (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
     b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
     b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
     b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
    a3 *
    (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
     b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
     b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
     b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

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
