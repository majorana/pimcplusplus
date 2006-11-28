#ifndef TRICUBIC_NUB_SPLINE_H
#define TRICUBIC_NUB_SPLINE_H

#include "BsplineHelper.h"
#include "Grid.h"

template<typename T, typename XGridType=LinearGrid, 
	 typename YGridType=XGridType, typename ZGridType=YGridType>
class TricubicNUBspline
{
private:
  NUBsplineBasis<XGridType> XBasis;
  NUBsplineBasis<YGridType> YBasis;
  NUBsplineBasis<ZGridType> ZBasis;
  int Nx, Ny, Nz;
  // The starting and ending values for the uniform grids
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  // The box dimensions and their inverses
  double Lx, LxInv, Ly, LyInv, Lz, LzInv;
  // The control points
  Array<T,3> P;

  TinyVector<double,3> Periodic;  
public:
  inline void Init (XGridType *xgrid, YGridType *ygrid, ZGridType *zgrid, Array<T,3> &data, 
		    BCType xbct=PERIODIC, BCType ybct=PERIODIC, BCType zbct=PERIODIC);
  inline T operator() (TinyVector<double,3> r);
  inline void Evaluate (TinyVector<double,3> r, double &val, 
			TinyVector<T,3> &grad);
  inline void Evaluate (TinyVector<double,3> r, T &val,
			TinyVector<T,3> &grad, T &laplacian);
  inline void Evaluate (TinyVector<double,3> r, T & val,
			TinyVector<T,3> &grad, 
			TinyMatrix<T,3,3> &secDerivs);
};


template<typename T, typename XGridType, typename YGridType, typename ZGridType>
void 
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::Init
(XGridType *xgrid, YGridType *ygrid, ZGridType *zgrid, Array<T,3> &data, 
 BCType xbc, BCType ybc, BCType zbc)
{
  // Set up 1D basis functions
  XBasis.Init (xgrid, xbc==PERIODIC);
  YBasis.Init (ygrid, ybc==PERIODIC);
  ZBasis.Init (zgrid, zbc==PERIODIC);
  Periodic[0] = xbc==PERIODIC ? 1.0 : 0.0;
  Periodic[1] = ybc==PERIODIC ? 1.0 : 0.0;
  Periodic[2] = zbc==PERIODIC ? 1.0 : 0.0;

  Nx = data.extent(0);
  Ny = data.extent(1);
  Nz = data.extent(2);

  int Mx, My, Mz;
  Mx = (xbc == PERIODIC) ? Nx+3 : Nx+2;
  My = (ybc == PERIODIC) ? Ny+3 : Ny+2;
  Mz = (zbc == PERIODIC) ? Nz+3 : Nz+2;
  assert (xgrid->NumPoints == Mx-2);
  assert (ygrid->NumPoints == My-2);
  assert (zgrid->NumPoints == Mz-2);

  P.resize(Mx, My, Mz);

  
  // Now solve interpolating equations
  
  TinyVector<T,4> lBC, rBC, dummy1, dummy2;
  ////////////////////
  // Do X direction //
  ////////////////////
  if (xbc == PERIODIC) {
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	SolvePeriodicInterp1D(XBasis, data(Range::all(), iy, iz), P(Range::all(),iy+1, iz+1));
  }
  else {
    if (xbc == FLAT) {
      XBasis (0   , dummy1, lBC);
      XBasis (Nx-1, dummy1, rBC);
    }
    else if (xbc == NATURAL) {
      XBasis (0   , dummy1, dummy2, lBC);
      XBasis (Nx-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	SolveDerivInterp1D(XBasis, data(Range::all(), iy, iz),  P(Range::all(), iy+1, iz+1), lBC, rBC);
  }
  ////////////////////
  // Do Y direction //
  ////////////////////
  if (ybc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++) 
	SolvePeriodicInterp1D(YBasis, P(ix, Range(1,Ny), iz), P(ix, Range::all(), iz));
  }
  else {
    if (ybc == FLAT) {
      YBasis (0,    dummy1, lBC);
      YBasis (Ny-1, dummy1, rBC);
    }
    else if (ybc == NATURAL) {
      YBasis (0,    dummy1, dummy2, lBC);
      YBasis (Ny-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iz=0; iz<Mz; iz++)
	SolveDerivInterp1D(YBasis, P(ix,Range(1,Ny), iz),  P(ix, Range::all(), iz), lBC, rBC);
  }
  ////////////////////
  // Do Z direction //
  ////////////////////
  if (zbc == PERIODIC) {
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++) 
	SolvePeriodicInterp1D(ZBasis, P(ix, iy, Range(1,Nz)), P(ix, iy, Range::all()));
  }
  else {
    if (zbc == FLAT) {
      ZBasis (0,    dummy1, lBC);
      ZBasis (Nz-1, dummy1, rBC);
    }
    else if (zbc == NATURAL) {
      ZBasis (0,    dummy1, dummy2, lBC);
      ZBasis (Nz-1, dummy1, dummy2, rBC);
    }
    lBC[3] = rBC[3] = 0.0;
    for (int ix=0; ix<Mx; ix++)
      for (int iy=0; iy<My; iy++)
	SolveDerivInterp1D(ZBasis, P(ix, iy, Range(1,Nz)),  P(ix, iy, Range::all()), lBC, rBC);
  }
}

template<typename T, typename XGridType, typename YGridType, typename ZGridType>
inline T 
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::operator() (TinyVector<double,3> r)
{
  TinyVector<double,4> a, b, c;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  // Return tensor product
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


#endif
