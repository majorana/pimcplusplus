#ifndef MULTI_TRICUBIC_SPLINE_H
#define MULTI_TRICUBIC_SPLINE_H

#include "Grid.h"
#include <cmath>
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

class MultiTricubicSpline
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
  void UpdateX (int source, int dest, int i);
  void UpdateY (int source, int dest, int i);
  void UpdateZ (int source, int dest, int i);
  bool UpToDate;
public:
  Array<TinyVector<double,8>,4> F;

  int Nx, Ny, Nz, N;
  Grid *Xgrid, *Ygrid, *Zgrid;
  TinyVector<Grid*,3> Grids;
  void Update();
  inline double operator()(int ix, int iy, int iz) const
  { return (F(ix,iy,ix)[0]); }
  inline double& operator()(int ix, int iy, int iz) 
  { UpToDate=false; return (F(ix,iy,iz)[0]); }
  inline void operator()(double x, double y, double z, Array<double,1> &vals);
  inline void d_dx      (double x, double y, double z, Array<double,1> &vals); 
  inline void d_dy      (double x, double y, double z, Array<double,1> &vals); 
  inline void d_dz      (double x, double y, double z, Array<double,1> &vals); 
  inline void d2_dx2    (double x, double y, double z, Array<double,1> &vals); 
  inline void d2_dy2    (double x, double y, double z, Array<double,1> &vals); 
  inline void d2_dz2    (double x, double y, double z, Array<double,1> &vals); 
  inline void d2_dxdy   (double x, double y, double z, Array<double,1> &vals); 
  inline void d2_dxdz   (double x, double y, double z, Array<double,1> &vals); 
  inline void d2_dydz   (double x, double y, double z, Array<double,1> &vals); 
  inline void Grad      (double x, double y, double z, Array<Vec3,  1> &grads);
  inline void Laplacian (double x, double y, double z, Array<double,1> &vals);

  MultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid, int n)
  {
    Xgrid = xgrid; Nx = xgrid->NumPoints;
    Ygrid = ygrid; Ny = ygrid->NumPoints;
    Zgrid = zgrid; Nz = zgrid->NumPoints;
    N = n;
    
    F.resize(Nx,Ny,Nz,N);
    UpToDate = false;
  }
  
  /// Copy constructor
  inline MultiTricubicSpline (const MultiTricubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline MultiTricubicSpline & operator= (MultiTricubicSpline &a);
  inline MultiTricubicSpline & operator= (MultiTricubicSpline a);


  inline void Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
		    const Array<double,4> &init);

  MultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid,
		   const Array<double,4> &init)
  {
    Init (xgrid, ygrid, zgrid, init);
  }
  MultiTricubicSpline() : UpToDate(false) 
  { /* Do nothing. */ }
};


inline MultiTricubicSpline::MultiTricubicSpline(const MultiTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
}


inline MultiTricubicSpline& MultiTricubicSpline::operator=(MultiTricubicSpline &a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}

inline MultiTricubicSpline& MultiTricubicSpline::operator=(MultiTricubicSpline a)
{
  F.resize(a.F.shape());
  F = a.F;
  Nx = a.Nx;
  Ny = a.Ny;
  Nz = a.Nz;
  N  = a.N;
  Xgrid = a.Xgrid;
  Ygrid = a.Ygrid;
  Zgrid = a.Zgrid;
  UpToDate = a.UpToDate;
  return (*this);
}



inline void MultiTricubicSpline::Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
				    const Array<double,4> &init)
{
  Xgrid = xgrid; Nx = xgrid->NumPoints;
  Ygrid = ygrid; Ny = ygrid->NumPoints;
  Zgrid = zgrid; Nz = zgrid->NumPoints;

  assert (init.extent(0) == Nx);
  assert (init.extent(1) == Ny);
  assert (init.extent(2) == Nz);
  N = init.extent(3);

  F.resize(Nx,Ny,Nz);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int i=0; i<N; i++)
	  F(ix,iy,iz,i)[0] = init(ix,iy,iz,i);
  UpToDate = false;
}




inline void MultiTricubicSpline::operator() (double x, double y, double z, 
					     Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

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

  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    vals(i) = 
      a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
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
  }
}


inline void 
MultiTricubicSpline::d_dx (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
  
    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    vals(i) = 
      da0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
  }
}



inline void
MultiTricubicSpline::d_dy (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);

  for (int i=0; i<N; i++) {
  
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    vals(i) = 
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
  }
}

inline void MultiTricubicSpline::d_dz (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

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

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  

  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    for (int i=0; i<N; i++)
      vals(i) = 
	a0*
	(b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
	 b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
	 b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
	 b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
	a1 *
	(b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
	 b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
	 b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
	 b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
	a2 *
	(b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
	 b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
	 b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
	 b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
	a3 *
	(b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
	 b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
	 b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
	 b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
  }
}

inline void
MultiTricubicSpline::Grad (double x, double y, double z, Array<Vec3,1> &grads)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,3> grad;
    grads(i)[0] = 
      da0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    grads(i)[1] = 
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    grads(i)[2] = 
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
  }
}


inline void 
MultiTricubicSpline::d2_dx2 (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double d2a0 = hinv*hinv*d2p1(u);
  double d2a1 = hinv*hinv*d2p2(u);
  double d2a2 = hinv*d2q1(u);
  double d2a3 = hinv*d2q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    double val = 
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

  }
}


inline void 
MultiTricubicSpline::d2_dy2 (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double d2b0 = kinv*kinv*d2p1(v);
  register double d2b1 = kinv*kinv*d2p2(v);
  register double d2b2 = kinv*d2q1(v);
  register double d2b3 = kinv*d2q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    vals(i) = 
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
  }
}


inline void 
MultiTricubicSpline::d2_dz2 (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

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

  register double d2c0 = linv*linv*d2p1(w);
  register double d2c1 = linv*linv*d2p2(w);
  register double d2c2 = linv*d2q1(w);
  register double d2c3 = linv*d2q2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    vals(i) = 
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
  }
}


inline void 
MultiTricubicSpline::Laplacian (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;

  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);
  double d2a0 = hinv*hinv*d2p1(u);
  double d2a1 = hinv*hinv*d2p2(u);
  double d2a2 = hinv*d2q1(u);
  double d2a3 = hinv*d2q2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);
  register double d2b0 = kinv*kinv*d2p1(v);
  register double d2b1 = kinv*kinv*d2p2(v);
  register double d2b2 = kinv*d2q1(v);
  register double d2b3 = kinv*d2q2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  register double d2c0 = linv*linv*d2p1(w);
  register double d2c1 = linv*linv*d2p2(w);
  register double d2c2 = linv*d2q1(w);
  register double d2c3 = linv*d2q2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    vals(i) = 
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
  
    vals(i) += 
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));

    vals(i) += 
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
  }
}

inline void 
MultiTricubicSpline::d2_dxdy (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double c0 = p1(w);
  register double c1 = p2(w);
  register double c2 = l*q1(w);
  register double c3 = l*q2(w);
  
  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    vals(i) = 
      da0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      da1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      da2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      da3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
  }
}

inline void 
MultiTricubicSpline::d2_dxdz (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double hinv = 1.0/h;
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double da0 = hinv*dp1(u);
  double da1 = hinv*dp2(u);
  double da2 = dq1(u);
  double da3 = dq2(u);

  register double b0 = p1(v);
  register double b1 = p2(v);
  register double b2 = k*q1(v);
  register double b3 = k*q2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);


  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    vals(i) = 
      da0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      da1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      da2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      da3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
  }
}


inline void 
MultiTricubicSpline::d2_dydz (double x, double y, double z, Array<double,1> &vals)
{
  if (!UpToDate)
    Update();

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);
  int iz = Zgrid->ReverseMap(z);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);
  iz = max(0,iz); iz = min(iz, Nz-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double kinv = 1.0/k;
  double l = (*Zgrid)(iz+1) - (*Zgrid)(iz);
  double linv = 1.0/l;
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  double w = (z - (*Zgrid)(iz))/l;

  double a0 = p1(u);
  double a1 = p2(u);
  double a2 = h*q1(u);
  double a3 = h*q2(u);

  register double db0 = kinv*dp1(v);
  register double db1 = kinv*dp2(v);
  register double db2 = dq1(v);
  register double db3 = dq2(v);

  register double dc0 = linv*dp1(w);
  register double dc1 = linv*dp2(w);
  register double dc2 = dq1(w);
  register double dc3 = dq2(w);

  for (int i=0; i<N; i++) {
    double& Y000 = F(ix,iy,iz,i)[0];      //   F
    double& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    double& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    double& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    double& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz

    double& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz

    double& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    vals(i) = 
      a0*
      (db0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       db1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       db2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       db3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (db0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       db1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       db2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       db3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (db0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       db1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       db2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       db3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (db0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       db1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       db2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       db3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
  }
}



#endif
