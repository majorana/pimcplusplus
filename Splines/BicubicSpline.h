#ifndef BICUBIC_SPLINE_H
#define BICUBIC_SPLINE_H

#include "CubicSpline.h"

/// This structure holds the required data for a grid point, which
/// includes the value of z at (xi,yj), dz/dx, dz/dy, and d^2/dxdy.
struct BCpoint
{
  double z;
  double dzdx;
  double dzdy;
  double d2zdxdy;
};

/// This class is used to interpolate a function z(x,y) in two
/// dimensions.  It is of cubic order in both directions.  It
/// currently uses only "natural" boundary conditions:  all second
/// derivatives are 0 at the boundary.
class BicubicSpline
{
private:
  /// Stores whether dzdx needs to be recalculated for each column
  Array<bool,1> XUpToDate;
  /// Stores whether dzdy needs to be recalculated for each row
  Array<bool,1> YUpToDate;
  /// Stores whether the d2zdxdy's need to be updated
  bool BiUpToDate;
  /// Update the iy column's dzdx's 
  void XUpdate(int iy);
  /// Update the ix row's dzdy's
  void YUpdate(int ix);
  /// Update all the d2zdxdy's
  void BiUpdate();
  /// Holds the BCpoint array containing z and its derivatives.
  Array<BCpoint,2> F;
public:
  int Nx, Ny;
  Grid *Xgrid, *Ygrid;

  /// Interpolates in y at a given row
  inline double  operator() (int ix,   double y);
  /// Interpolates in x at a given column
  inline double  operator() (double x, int iy);
  /// Bicubic interpolation in x and y
  inline double  operator() (double x, double y);
  /// Returns z(x_{ix}, y_{iy})
  inline double  operator() (int ix, int iy) const;
  /// Returns a reference to z(x_{ix}, y_{iy}) that can be used as an
  /// L-value.
  inline double& operator() (int ix, int iy);
  inline double  Deriv      (int ix,   double y);
  inline double  Deriv      (double x, int iy);
  inline double  xDeriv     (int ix, int iy);
  inline double  yDeriv     (int ix, int iy);
  inline double  Deriv2     (int ix,   double y);
  inline double  Deriv2     (double x, int iy);
  inline double  Deriv3     (int ix,   double y);
  inline double  Deriv3     (double x, int iy);
  /// Initialize the bicubic spline with the given grids and data
  inline void Init (Grid *xgrid, Grid *ygrid, Array<double,2> &f);
};

inline void BicubicSpline::Init(Grid *xgrid, Grid *ygrid, Array<double,2> &f)
{
  Nx = xgrid->NumPoints;
  Ny = ygrid->NumPoints;
  Xgrid = xgrid; Ygrid = ygrid;
  assert (f.rows() == Nx);
  assert (f.cols() == Ny);
  XUpToDate.resize(Ny);
  YUpToDate.resize(Nx);
  F.resize(f.rows(),f.cols());
  for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
      F(i,j).z = f(i,j);
  for (int i=0; i<Ny; i++)
    XUpdate(i);
  for (int i=0; i<Nx; i++)
    YUpdate(i);

  BiUpdate();
}

inline double BicubicSpline::operator() (int ix, int iy) const
{ return (F(ix,iy).z);}

inline double& BicubicSpline::operator() (int ix, int iy)
{ 
  XUpToDate(iy) = false;
  YUpToDate(ix) = false;
  BiUpToDate = false;
  return (F(ix,iy).z);
}

inline double BicubicSpline::operator() (double x, int iy)
{
  if (!XUpToDate(iy))
    XUpdate(iy);
  
  int ix = Xgrid->ReverseMap(x);
  ix = max(0,ix);
  ix = min(ix, Nx-2);
  
  double t = (x - (*Xgrid)(ix))/((*Xgrid)(ix+1) - (*Xgrid)(ix));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);

  return (F(ix,iy).z*p1 + F(ix+1,iy).z*p2 + 
	  h*(F(ix,iy).dzdx*q1 + F(ix+1,iy).dzdx*q2));
}


inline double BicubicSpline::operator() (int ix, double y)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);

  return (F(ix,iy).z*p1 + F(ix,iy+1).z*p2 + 
	  h*(F(ix,iy).dzdy*q1 + F(ix,iy+1).dzdy*q2));
}


inline double BicubicSpline::Deriv (int ix, double y)
{
  cerr << "BicubicSpline::Deriv Not implemented yet!\n";
  abort();
}

inline double BicubicSpline::Deriv2 (int ix, double y)
{
  cerr << "BicubicSpline::Deriv2 Not implemented yet!\n";
  abort();
}

inline double BicubicSpline::Deriv3 (int ix, double y)
{
  cerr << "BicubicSpline::Deriv3 Not implemented yet!\n";
  abort();
}


inline double p1(double t)
{ return ((t-1.0)*(t-1.0)*(1.0+2.0*t)); }

inline double p2(double t)
{ return (t*t*(3.0-2.0*t)); }

inline double q1(double t)
{ return (t*(t-1.0)*(t-1.0)); }

inline double q2(double t)
{ return (t*t*(t-1.0)); }

inline double BicubicSpline::operator() (double x, double y)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = F(ix,iy).z;
  Z(0,1) = F(ix,iy+1).z;
  Z(0,2) = F(ix,iy).dzdy;
  Z(0,3) = F(ix,iy+1).dzdy;
  Z(1,0) = F(ix+1,iy).z;
  Z(1,1) = F(ix+1,iy+1).z;
  Z(1,2) = F(ix+1,iy).dzdy;
  Z(1,3) = F(ix+1,iy+1).dzdy;
  Z(2,0) = F(ix,iy).dzdx;
  Z(2,1) = F(ix,iy+1).dzdx;
  Z(2,2) = F(ix,iy).d2zdxdy;
  Z(2,3) = F(ix,iy+1).d2zdxdy;
  Z(3,0) = F(ix+1,iy).dzdx;
  Z(3,1) = F(ix+1,iy+1).dzdx;
  Z(3,2) = F(ix+1,iy).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}



class MultiBicubicSpline
{
private:
  Array<bool,1> XUpToDate;
  Array<bool,1> YUpToDate;
  bool BiUpToDate;
  void XUpdate(int ix);
  void YUpdate(int iy);
  void BiUpdate();
  Array<BCpoint,3> F;
public:
  int Nx, Ny, Nz;
  Grid *Xgrid, *Ygrid;

  inline double  operator() (int ix,   double y, int iz);
  inline double  operator() (double x, int iy, int iz);
  inline double  operator() (double x, double y, int iz);
  inline double  operator() (int ix, int iy, int iz) const;
  inline double& operator() (int ix, int iy, int iz);
  inline void    operator() (double x, double y, Array<double,1> &z);
  inline double  Deriv      (int ix,   double y, int iz);
  inline double  Deriv      (double x, int iy, int iz);
  inline double  xDeriv     (int ix, int iy, int iz);
  inline double  yDeriv     (int ix, int iy, int iz);
  inline double  Deriv2     (int ix,   double y, int iz);
  inline double  Deriv2     (double x, int iy, int iz);
  inline double  Deriv3     (int ix,   double y, int iz);
  inline double  Deriv3     (double x, int iy, int iz);
  inline void Init (Grid *xgrid, Grid *ygrid, Array<double,3> &f);
};

inline void 
MultiBicubicSpline::Init(Grid *xgrid, Grid *ygrid, Array<double,3> &f)
{
  Nx = xgrid->NumPoints;
  Ny = ygrid->NumPoints;
  Nz = f.extent(2);
  Xgrid = xgrid; Ygrid = ygrid;
  assert (f.rows() == Nx);
  assert (f.cols() == Ny);
  XUpToDate.resize(Ny);
  YUpToDate.resize(Nx);
  F.resize(Nx, Ny, Nz);
  for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
      for (int k=0; k<Nz; k++)
	F(i,j,k).z = f(i,j,k);
  for (int i=0; i<Ny; i++)
    XUpdate(i);
  for (int i=0; i<Nx; i++)
    YUpdate(i);

  BiUpdate();
}


inline double MultiBicubicSpline::operator() (int ix, int iy, int iz) const
{ return (F(ix,iy,iz).z);}

inline double& MultiBicubicSpline::operator() (int ix, int iy, int iz)
{ 
  XUpToDate(iy) = false;
  YUpToDate(ix) = false;
  BiUpToDate = false;
  return (F(ix,iy,iz).z);
}


inline double MultiBicubicSpline::operator() (double x, int iy, int iz)
{
  if (!XUpToDate(iy))
    XUpdate(iy);
  
  int ix = Xgrid->ReverseMap(x);
  ix = max(0,ix);
  ix = min(ix, Nx-2);
  
  double t = (x - (*Xgrid)(ix))/((*Xgrid)(ix+1) - (*Xgrid)(ix));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);

  return (F(ix,iy,iz).z*p1 + F(ix+1,iy,iz).z*p2 + 
	  h*(F(ix,iy,iz).dzdx*q1 + F(ix+1,iy,iz).dzdx*q2));
}


inline double MultiBicubicSpline::operator() (int ix, double y, int iz)
{
  if (!YUpToDate(ix))
    YUpdate(ix);
  
  int iy = Ygrid->ReverseMap(y);
  iy = max(0,iy);
  iy = min(iy, Ny-2);
  
  double t = (y - (*Ygrid)(iy))/((*Ygrid)(iy+1) - (*Ygrid)(iy));
  double tm1 = t - 1.0;
  double p1 = tm1*tm1*(1.0+2.0*t);
  double q1 = t*tm1*tm1;
  double p2 = t*t*(3.0-2.0*t);
  double q2 = t*t*tm1;
  double h = (*Ygrid)(iy+1) - (*Ygrid)(iy);

  return (F(iy,ix,iz).z*p1 + F(iy+1,ix,iz).z*p2 + 
	  h*(F(iy,ix,iz).dzdy*q1 + F(iy+1,ix,iz).dzdy*q2));
}

inline double MultiBicubicSpline::operator() (double x, double y, int iz)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  Z(0,0) = F(ix,iy,iz).z;
  Z(0,1) = F(ix,iy+1,iz).z;
  Z(0,2) = F(ix,iy,iz).dzdy;
  Z(0,3) = F(ix,iy+1,iz).dzdy;
  Z(1,0) = F(ix+1,iy,iz).z;
  Z(1,1) = F(ix+1,iy+1,iz).z;
  Z(1,2) = F(ix+1,iy,iz).dzdy;
  Z(1,3) = F(ix+1,iy+1,iz).dzdy;
  Z(2,0) = F(ix,iy,iz).dzdx;
  Z(2,1) = F(ix,iy+1,iz).dzdx;
  Z(2,2) = F(ix,iy,iz).d2zdxdy;
  Z(2,3) = F(ix,iy+1,iz).d2zdxdy;
  Z(3,0) = F(ix+1,iy,iz).dzdx;
  Z(3,1) = F(ix+1,iy+1,iz).dzdx;
  Z(3,2) = F(ix+1,iy,iz).d2zdxdy;
  Z(3,3) = F(ix+1,iy+1,iz).d2zdxdy;
  
  double val = 0.0;
  for (int m=0; m<4; m++) {
    double Zb_m = 0.0;
    for (int n=0; n<4; n++)
      Zb_m += Z(m,n) * b(n);
    val += Zb_m * a(m);
  }
  return (val);
}



inline void MultiBicubicSpline::operator() (double x, double y, 
					    Array<double,1> &z)
{
  if (!BiUpToDate)
    BiUpdate();
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Xgrid->ReverseMap(x);  
  int iy = Ygrid->ReverseMap(y);

  ix = max(0,ix); ix = min(ix, Nx-2);
  iy = max(0,iy); iy = min(iy, Ny-2);

  double h = (*Xgrid)(ix+1) - (*Xgrid)(ix);
  double k = (*Ygrid)(iy+1) - (*Ygrid)(iy);
  double u = (x - (*Xgrid)(ix))/h;
  double v = (y - (*Ygrid)(iy))/k;
  a(0) = p1(u);
  a(1) = p2(u);
  a(2) = h*q1(u);
  a(3) = h*q2(u);

  b(0) = p1(v);
  b(1) = p2(v);
  b(2) = k*q1(v);
  b(3) = k*q2(v);
  
  for (int iz=0; iz<Nz; iz++) {
    Z(0,0) = F(ix,iy,iz).z;
    Z(0,1) = F(ix,iy+1,iz).z;
    Z(0,2) = F(ix,iy,iz).dzdy;
    Z(0,3) = F(ix,iy+1,iz).dzdy;
    Z(1,0) = F(ix+1,iy,iz).z;
    Z(1,1) = F(ix+1,iy+1,iz).z;
    Z(1,2) = F(ix+1,iy,iz).dzdy;
    Z(1,3) = F(ix+1,iy+1,iz).dzdy;
    Z(2,0) = F(ix,iy,iz).dzdx;
    Z(2,1) = F(ix,iy+1,iz).dzdx;
    Z(2,2) = F(ix,iy,iz).d2zdxdy;
    Z(2,3) = F(ix,iy+1,iz).d2zdxdy;
    Z(3,0) = F(ix+1,iy,iz).dzdx;
    Z(3,1) = F(ix+1,iy+1,iz).dzdx;
    Z(3,2) = F(ix+1,iy,iz).d2zdxdy;
    Z(3,3) = F(ix+1,iy+1,iz).d2zdxdy;
  
    z(iz) = 0.0;
    for (int m=0; m<4; m++) {
      double Zb_m = 0.0;
      for (int n=0; n<4; n++)
	Zb_m += Z(m,n) * b(n);
      z(iz) += Zb_m * a(m);
    }
  }
}



#endif
