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


#endif
