#ifndef BICUBIC_SPLINE_H
#define BICUBIC_SPLINE_H

#include "CubicSpline.h"

class BicubicSpline
{
private:
  Array<bool,1> XUpToDate;
  Array<bool,1> YUpToDate;
  bool BiUpToDate;
  Array<double,2> F;
  Array<double,2> d2Fdx2;
  Array<double,2> d2Fdy2;
  void XUpdate(int ix);
  void YUpdate(int iy);
  void BiUpdate();
  Array<double,1> XStartDeriv, XEndDeriv, YStartDeriv, YEndDeriv;
public:
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

inline double BicubicSpline::operator() (int ix, int iy) const
{ return (F(ix,iy));}

inline double& BicubicSpline::operator() (int ix, int iy)
{ 
  XUpToDate(iy) = false;
  YUpToDate(ix) = false;
  return (F(ix,iy));
}



#endif
