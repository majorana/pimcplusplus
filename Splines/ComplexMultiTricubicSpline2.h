#ifndef COMPLEX_MULTI_TRICUBIC_SPLINE_H
#define COMPLEX_MULTI_TRICUBIC_SPLINE_H

#include "Grid.h"
#include <cmath>
//#include <blitz/array.h>
//using namespace blitz;

/// Each point of F contains:
/// 0)  real [F(x,y,z)]
/// 1)  real [dF/dx]
/// 2)  real [dF/dy]
/// 3)  real [dF/dz]
/// 4)  real [d2F/dxdy]
/// 5)  real [d2F/dxdz]
/// 6)  real [d2F/dydz]
/// 7)  real [d3F/dxdydz]
/// 8)  imag [F(x,y,z)]
/// 9)  imag [dF/dx]
/// 10) imag [dF/dy]
/// 11) imag [dF/dz]
/// 12) imag [d2F/dxdy]
/// 13) imag [d2F/dxdz]
/// 14) imag [d2F/dydz]
/// 15) imag [d3F/dxdydz]


class ComplexMultiTricubicSpline
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
  void UpdateXPeriodic (int source, int dest, int i);
  void UpdateYPeriodic (int source, int dest, int i);
  void UpdateZPeriodic (int source, int dest, int i);
  bool UpToDate, Periodic;
public:
  Array<TinyVector<TinyVector<double,2>,8>,4> F;

  int Nx, Ny, Nz, N;
  Grid *Xgrid, *Ygrid, *Zgrid;
  TinyVector<Grid*,3> Grids;
  void Update();
  inline complex<double> operator()(int ix, int iy, int iz, int i) const
  { return complex<double> (F(ix,iy,ix, i)[0][0], F(ix,iy,iz,i)[0][1]); }
  inline void Set (int ix, int iy, int iz, int i, complex<double> val)
    { UpToDate=false; F(ix, iy, iz, i)[0][0] = val.real(); F(ix,iy,iz,i)[0][1] = val.imag(); }
  inline void operator()(double x, double y, double z, 
			 Array<complex<double>,1> &vals);
  inline void d_dx      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d_dy      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d_dz      (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dx2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dy2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dz2    (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dxdy   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dxdz   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void d2_dydz   (double x, double y, double z, 
			 Array<complex<double>,1> &vals); 
  inline void Grad      (double x, double y, double z, 
			 Array<cVec3,1> &grads);
  inline void ValGrad   (double x, double y, double z, 
			 Array<complex<double>,1> &vals, 
			 Array<cVec3,  1> &grads);
  inline void Laplacian (double x, double y, double z, 
			 Array<complex<double>,1> &vals);
  inline void* Data () { return F.data(); }

  ComplexMultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid, int n)
  {
    Xgrid = xgrid; Nx = xgrid->NumPoints;
    Ygrid = ygrid; Ny = ygrid->NumPoints;
    Zgrid = zgrid; Nz = zgrid->NumPoints;
    N = n;
    
    F.resize(Nx,Ny,Nz,N);
    UpToDate = false;
  }
  
  /// Copy constructor
  inline ComplexMultiTricubicSpline (const ComplexMultiTricubicSpline &a);

  /// Assigment operator -- necessary for array resizeAndPreserve
  inline ComplexMultiTricubicSpline & operator= (ComplexMultiTricubicSpline &a);
  inline ComplexMultiTricubicSpline & operator= (ComplexMultiTricubicSpline a);


  inline void Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
		    const Array<complex<double>,4> &init, bool periodic=false);

  ComplexMultiTricubicSpline(Grid *xgrid, Grid *ygrid, Grid *zgrid,
		      const Array<complex<double>,4> &init, bool periodic=false)
  {
    Init (xgrid, ygrid, zgrid, init, periodic);
  }
  ComplexMultiTricubicSpline() : UpToDate(false) 
  { /* Do nothing. */ }
};

inline 
ComplexMultiTricubicSpline::ComplexMultiTricubicSpline
(const ComplexMultiTricubicSpline &a)
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


inline 
ComplexMultiTricubicSpline& ComplexMultiTricubicSpline::operator=
(ComplexMultiTricubicSpline &a)
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

inline ComplexMultiTricubicSpline& 
ComplexMultiTricubicSpline::operator=(ComplexMultiTricubicSpline a)
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



inline void 
ComplexMultiTricubicSpline::Init (Grid *xgrid, Grid *ygrid, Grid *zgrid,
				  const Array<complex<double>,4> &init,
				  bool periodic)
{
  Periodic = periodic;

  Xgrid = xgrid; Nx = xgrid->NumPoints;
  Ygrid = ygrid; Ny = ygrid->NumPoints;
  Zgrid = zgrid; Nz = zgrid->NumPoints;

  assert (init.extent(0) == Nx);
  assert (init.extent(1) == Ny);
  assert (init.extent(2) == Nz);
  N = init.extent(3);

  F.resize(Nx,Ny,Nz,N);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int i=0; i<N; i++) {
	  F(ix,iy,iz,i)[0][0] = init(ix,iy,iz,i).real();
	  F(ix,iy,iz,i)[8][1] = init(ix,iy,iz,i).imag();
	}
  UpToDate = false;
  Update();
}




inline void ComplexMultiTricubicSpline::operator()
  (double x, double y, double z, Array<complex<double>,1> &vals)
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
    //////////////////
    /// Real parts ///
    //////////////////
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2> Val;
    Val = 
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
    vals(i) = complex<double>(Val[0],Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d_dx (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    
    TinyVector<double,2> Val;
    Val = 
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

    vals(i) = complex<double>(Val[0],Val[1]);
  }
}



inline void
ComplexMultiTricubicSpline::d_dy (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val =
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

    vals(i) = complex<double> (Val[0], Val[1]);
  }
}

inline void 
ComplexMultiTricubicSpline::d_dz (double x, double y, double z, 
				  Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val = 
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

    vals(i) = complex<double> (Val[0], Val[1]);
  }
}

inline void
ComplexMultiTricubicSpline::Grad (double x, double y, double z, 
				  Array<cVec3,1> &grads)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<TinyVector<double,2>,3> Grad;
    Grad[0] = 
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

    Grad[1] = 
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

    Grad[2] = 
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

    grads(i)[0] = complex<double> (Grad[0][0], Grad[0][1]);
    grads(i)[1] = complex<double> (Grad[1][0], Grad[1][1]);
    grads(i)[2] = complex<double> (Grad[2][0], Grad[2][1]);
  }
}

/// Returns the value and computes the gradient
inline void
ComplexMultiTricubicSpline::ValGrad(double x, double y, double z, 
				    Array<complex<double>,1> &vals, 
				    Array<cVec3,1> &grads)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,2> Val;
    TinyVector<TinyVector<double,2>,3> Grad;
    Val = 
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

    vals(i) = complex<double> (Val[0], Val[1]);


    Grad[0] = 
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

    Grad[1] = 
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

    Grad[2] = 
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
    grads(i)[0] = complex<double>(Grad[0][0], Grad[0][1]);
    grads(i)[1] = complex<double>(Grad[1][0], Grad[1][1]);
    grads(i)[2] = complex<double>(Grad[2][0], Grad[2][1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dx2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,2> Val =
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

    vals(i) = complex<double> (Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dy2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz


    TinyVector<double,2> Val = 
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

      vals(i) = complex<double> (Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dz2 (double x, double y, double z, 
				    Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
  
    TinyVector<double,2> Val = 
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

    vals(i) = complex<double>(Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::Laplacian (double x, double y, double z, 
				       Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val =
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
  
    Val +=
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

    Val += 
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

    vals(i) = complex<double>(Val[0],Val[1]);
  }
}

inline void 
ComplexMultiTricubicSpline::d2_dxdy (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2> Val =
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

    vals(i) = complex<double>(Val[0], Val[1]);
  }
}

inline void 
ComplexMultiTricubicSpline::d2_dxdz (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val = 
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
    vals(i) = complex<double>(Val[0], Val[1]);
  }
}


inline void 
ComplexMultiTricubicSpline::d2_dydz (double x, double y, double z, 
				     Array<complex<double>,1> &vals)
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
    TinyVector<double,2>& Y000 = F(ix,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y001 = F(ix,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y002 = F(ix,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y003 = F(ix,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y010 = F(ix,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y011 = F(ix,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y012 = F(ix,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y013 = F(ix,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y020 = F(ix,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y021 = F(ix,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y022 = F(ix,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y023 = F(ix,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y030 = F(ix,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y031 = F(ix,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y032 = F(ix,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y033 = F(ix,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y100 = F(ix+1,iy,iz,i)[0];      //   F
    TinyVector<double,2>& Y101 = F(ix+1,iy,iz+1,i)[0];    //   F
    TinyVector<double,2>& Y102 = F(ix+1,iy,iz,i)[3];      //  dF/dz
    TinyVector<double,2>& Y103 = F(ix+1,iy,iz+1,i)[3];    //  dF/dz
    TinyVector<double,2>& Y110 = F(ix+1,iy+1,iz,i)[0];    //   F
    TinyVector<double,2>& Y111 = F(ix+1,iy+1,iz+1,i)[0];  //   F
    TinyVector<double,2>& Y112 = F(ix+1,iy+1,iz,i)[3];    //  dF/dz
    TinyVector<double,2>& Y113 = F(ix+1,iy+1,iz+1,i)[3];  //  dF/dz
    TinyVector<double,2>& Y120 = F(ix+1,iy,iz,i)[2];      //  dF/dy
    TinyVector<double,2>& Y121 = F(ix+1,iy,iz+1,i)[2];    //  dF/dy
    TinyVector<double,2>& Y122 = F(ix+1,iy,iz,i)[6];      // d2F/dydz
    TinyVector<double,2>& Y123 = F(ix+1,iy,iz+1,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y130 = F(ix+1,iy+1,iz,i)[2];    //  dF/dy
    TinyVector<double,2>& Y131 = F(ix+1,iy+1,iz+1,i)[2];  //  dF/dy
    TinyVector<double,2>& Y132 = F(ix+1,iy+1,iz,i)[6];    // d2F/dydz
    TinyVector<double,2>& Y133 = F(ix+1,iy+1,iz+1,i)[6];  // d2F/dydz
    
    TinyVector<double,2>& Y200 = F(ix,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y201 = F(ix,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y202 = F(ix,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y203 = F(ix,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y210 = F(ix,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y211 = F(ix,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y212 = F(ix,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y213 = F(ix,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y220 = F(ix,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y221 = F(ix,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y222 = F(ix,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y223 = F(ix,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y230 = F(ix,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y231 = F(ix,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y232 = F(ix,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y233 = F(ix,iy+1,iz+1,i)[7];  // d3F/dxdydz
    
    TinyVector<double,2>& Y300 = F(ix+1,iy,iz,i)[1];      //  dF/dx
    TinyVector<double,2>& Y301 = F(ix+1,iy,iz+1,i)[1];    //  dF/dx
    TinyVector<double,2>& Y302 = F(ix+1,iy,iz,i)[5];      // d2F/dxdz
    TinyVector<double,2>& Y303 = F(ix+1,iy,iz+1,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y310 = F(ix+1,iy+1,iz,i)[1];    //  dF/dx
    TinyVector<double,2>& Y311 = F(ix+1,iy+1,iz+1,i)[1];  //  dF/dx
    TinyVector<double,2>& Y312 = F(ix+1,iy+1,iz,i)[5];    // d2F/dxdz
    TinyVector<double,2>& Y313 = F(ix+1,iy+1,iz+1,i)[5];  // d2F/dxdz
    TinyVector<double,2>& Y320 = F(ix+1,iy,iz,i)[4];      // d2F/dxdy
    TinyVector<double,2>& Y321 = F(ix+1,iy,iz+1,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y322 = F(ix+1,iy,iz,i)[7];      // d3F/dxdydz
    TinyVector<double,2>& Y323 = F(ix+1,iy,iz+1,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y330 = F(ix+1,iy+1,iz,i)[4];    // d2F/dxdy
    TinyVector<double,2>& Y331 = F(ix+1,iy+1,iz+1,i)[4];  // d2F/dxdy
    TinyVector<double,2>& Y332 = F(ix+1,iy+1,iz,i)[7];    // d3F/dxdydz
    TinyVector<double,2>& Y333 = F(ix+1,iy+1,iz+1,i)[7];  // d3F/dxdydz

    TinyVector<double,2> Val = 
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
    vals(i) = complex<double>(Val[0], Val[1]);
  }
}

/// This replicates the first points of the 3D array to the last
/// points, making the function periodic
void MakePeriodic(Array<complex<double>,4> &A);

#endif
