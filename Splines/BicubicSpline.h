#ifndef BICUBIC_SPLINE_H
#define BICUBIC_SPLINE_H

#include "CubicSpline.h"

struct BCpoint
{
  double z;
  double dzdx;
  double dzdy;
  double d2zdxdy;
};


class BicubicSpline
{
private:
  Array<bool,1> XUpToDate;
  Array<bool,1> YUpToDate;
  bool BiUpToDate;
  void XUpdate(int ix);
  void YUpdate(int iy);
  void BiUpdate();
public:
  Array<BCpoint,2> F;
  int Nx, Ny;
  Grid *Xgrid, *Ygrid;

  inline double  operator() (int ix,   double y);
  inline double  operator() (double x, int iy);
  inline double  operator() (double x, double y);
  inline double  operator() (int ix, int iy) const;
  inline double& operator() (int ix, int iy);
  inline double  Deriv      (int ix,   double y);
  inline double  Deriv      (double x, int iy);
  inline double  xDeriv     (int ix, int iy);
  inline double  yDeriv     (int ix, int iy);
  inline double  Deriv2     (int ix,   double y);
  inline double  Deriv2     (double x, int iy);
  inline double  Deriv3     (int ix,   double y);
  inline double  Deriv3     (double x, int iy);
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

  return (F(iy,ix).z*p1 + F(iy+1,ix).z*p2 + 
	  h*(F(iy,ix).dzdy*q1 + F(iy+1,ix).dzdy*q2));
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
  TinyMatrix<double,4,4> Z;
  TinyVector<double,4> a, b;

  int ix = Ygrid->ReverseMap(x);  
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

#endif
