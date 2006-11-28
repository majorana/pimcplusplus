#ifndef TRICUBIC_NUB_SPLINE_H
#define TRICUBIC_NUB_SPLINE_H

#include "BsplineHelper.h"
#include "NUBsplineBasis.h"
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
  inline T operator() (TinyVector<double,3> r) const;
  inline void Evaluate (TinyVector<double,3> r, T &val, 
			TinyVector<T,3> &grad) const;
  inline void Evaluate (TinyVector<double,3> r, T &val,
			TinyVector<T,3> &grad, T &laplacian) const;
  inline void Evaluate (TinyVector<double,3> r, T & val,
			TinyVector<T,3> &grad, 
			TinyMatrix<T,3,3> &secDerivs) const;
};


void Duplicate (TinyVector<double,4> source,
		TinyVector<complex<double>,4> dest)
{
  dest[0] = complex<double>(source[0], source[0]);
  dest[1] = complex<double>(source[1], source[1]);
  dest[2] = complex<double>(source[2], source[2]);
  dest[3] = complex<double>(source[3], source[3]);
}

void Duplicate (TinyVector<double,4> source,
		TinyVector<double,4> dest)
{
  dest = source;
}



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
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::operator() (TinyVector<double,3> r) const
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

template<typename T, typename XGridType, typename YGridType, typename ZGridType>
inline void
TricubicNUBspline<T,XGridType,YGridType,ZGridType>::Evaluate (TinyVector<double,3> r,
							      T &val,
							      TinyVector<T,3> &grad) const
{
  TinyVector<double,4> a, da, b, db, c, dc;
  // Evaluate 1D basis functions
  int ix0 = XBasis(r[0], a, da); int ix1=ix0+1; int ix2=ix0+2; int ix3=ix0+3;
  int iy0 = YBasis(r[1], b, db); int iy1=iy0+1; int iy2=iy0+2; int iy3=iy0+3;
  int iz0 = ZBasis(r[2], c, dc); int iz1=iz0+1; int iz2=iz0+2; int iz3=iz0+3;

  T Pi[64], cP[16], bcP[4];
  Pi[ 0]=P(ix0,iy0,iz0);  Pi[ 1]=P(ix0,iy0,iz1);  Pi[ 2]=P(ix0,iy0,iz2);  Pi[ 3]=P(ix0,iy0,iz3);
  Pi[ 4]=P(ix0,iy1,iz0);  Pi[ 5]=P(ix0,iy1,iz1);  Pi[ 6]=P(ix0,iy1,iz2);  Pi[ 7]=P(ix0,iy1,iz3);
  Pi[ 8]=P(ix0,iy2,iz0);  Pi[ 9]=P(ix0,iy2,iz1);  Pi[10]=P(ix0,iy2,iz2);  Pi[11]=P(ix0,iy2,iz3);
  Pi[12]=P(ix0,iy3,iz0);  Pi[13]=P(ix0,iy3,iz1);  Pi[14]=P(ix0,iy3,iz2);  Pi[15]=P(ix0,iy3,iz3);
  Pi[16]=P(ix1,iy0,iz0);  Pi[17]=P(ix1,iy0,iz1);  Pi[18]=P(ix1,iy0,iz2);  Pi[19]=P(ix1,iy0,iz3);
  Pi[20]=P(ix1,iy1,iz0);  Pi[21]=P(ix1,iy1,iz1);  Pi[22]=P(ix1,iy1,iz2);  Pi[23]=P(ix1,iy1,iz3);
  Pi[24]=P(ix1,iy2,iz0);  Pi[25]=P(ix1,iy2,iz1);  Pi[26]=P(ix1,iy2,iz2);  Pi[27]=P(ix1,iy2,iz3);
  Pi[28]=P(ix1,iy3,iz0);  Pi[29]=P(ix1,iy3,iz1);  Pi[30]=P(ix1,iy3,iz2);  Pi[31]=P(ix1,iy3,iz3);
  Pi[32]=P(ix2,iy0,iz0);  Pi[33]=P(ix2,iy0,iz1);  Pi[34]=P(ix2,iy0,iz2);  Pi[35]=P(ix2,iy0,iz3);
  Pi[36]=P(ix2,iy1,iz0);  Pi[37]=P(ix2,iy1,iz1);  Pi[38]=P(ix2,iy1,iz2);  Pi[39]=P(ix2,iy1,iz3);
  Pi[40]=P(ix2,iy2,iz0);  Pi[41]=P(ix2,iy2,iz1);  Pi[42]=P(ix2,iy2,iz2);  Pi[43]=P(ix2,iy2,iz3);
  Pi[44]=P(ix2,iy3,iz0);  Pi[45]=P(ix2,iy3,iz1);  Pi[46]=P(ix2,iy3,iz2);  Pi[47]=P(ix2,iy3,iz3);
  Pi[48]=P(ix3,iy0,iz0);  Pi[49]=P(ix3,iy0,iz1);  Pi[50]=P(ix3,iy0,iz2);  Pi[51]=P(ix3,iy0,iz3);
  Pi[52]=P(ix3,iy1,iz0);  Pi[53]=P(ix3,iy1,iz1);  Pi[54]=P(ix3,iy1,iz2);  Pi[55]=P(ix3,iy1,iz3);
  Pi[56]=P(ix3,iy2,iz0);  Pi[57]=P(ix3,iy2,iz1);  Pi[58]=P(ix3,iy2,iz2);  Pi[59]=P(ix3,iy2,iz3);
  Pi[60]=P(ix3,iy3,iz0);  Pi[61]=P(ix3,iy3,iz1);  Pi[62]=P(ix3,iy3,iz2);  Pi[63]=P(ix3,iy3,iz3);

  cP[ 0] = c[0]*Pi[ 0] + c[1]*Pi[ 1] + c[2]*Pi[ 2] + c[3]*Pi[ 3];
  cP[ 1] = c[0]*Pi[ 4] + c[1]*Pi[ 5] + c[2]*Pi[ 6] + c[3]*Pi[ 7];
  cP[ 2] = c[0]*Pi[ 8] + c[1]*Pi[ 9] + c[2]*Pi[10] + c[3]*Pi[11];
  cP[ 3] = c[0]*Pi[12] + c[1]*Pi[13] + c[2]*Pi[14] + c[3]*Pi[15];
  cP[ 4] = c[0]*Pi[16] + c[1]*Pi[17] + c[2]*Pi[18] + c[3]*Pi[19];
  cP[ 5] = c[0]*Pi[20] + c[1]*Pi[21] + c[2]*Pi[22] + c[3]*Pi[23];
  cP[ 6] = c[0]*Pi[24] + c[1]*Pi[25] + c[2]*Pi[26] + c[3]*Pi[27];
  cP[ 7] = c[0]*Pi[28] + c[1]*Pi[29] + c[2]*Pi[30] + c[3]*Pi[31];
  cP[ 8] = c[0]*Pi[32] + c[1]*Pi[33] + c[2]*Pi[34] + c[3]*Pi[35];
  cP[ 9] = c[0]*Pi[36] + c[1]*Pi[37] + c[2]*Pi[38] + c[3]*Pi[39];
  cP[10] = c[0]*Pi[40] + c[1]*Pi[41] + c[2]*Pi[42] + c[3]*Pi[43];
  cP[11] = c[0]*Pi[44] + c[1]*Pi[45] + c[2]*Pi[46] + c[3]*Pi[47];
  cP[12] = c[0]*Pi[48] + c[1]*Pi[49] + c[2]*Pi[50] + c[3]*Pi[51];
  cP[13] = c[0]*Pi[52] + c[1]*Pi[53] + c[2]*Pi[54] + c[3]*Pi[55];
  cP[14] = c[0]*Pi[56] + c[1]*Pi[57] + c[2]*Pi[58] + c[3]*Pi[59];
  cP[15] = c[0]*Pi[60] + c[1]*Pi[61] + c[2]*Pi[62] + c[3]*Pi[63];

  bcP[0] = b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3];
  bcP[1] = b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7];
  bcP[2] = b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11];
  bcP[3] = b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15];

  val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];
  grad[0] = da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3];
  grad[1] = (a[0]*(db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]) +
	     a[1]*(db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]) +
	     a[2]*(db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]) +
	     a[3]*(db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]));
  grad[2] = 
    a[0]*(b[0]*(dc[0]*Pi[ 0] + dc[1]*Pi[ 1] + dc[2]*Pi[ 2] + dc[3]*Pi[ 3])+
	  b[1]*(dc[0]*Pi[ 4] + dc[1]*Pi[ 5] + dc[2]*Pi[ 6] + dc[3]*Pi[ 7])+
	  b[2]*(dc[0]*Pi[ 8] + dc[1]*Pi[ 9] + dc[2]*Pi[10] + dc[3]*Pi[11])+
	  b[3]*(dc[0]*Pi[12] + dc[1]*Pi[13] + dc[2]*Pi[14] + dc[3]*Pi[15]))+
    a[1]*(b[0]*(dc[0]*Pi[16] + dc[1]*Pi[17] + dc[2]*Pi[18] + dc[3]*Pi[19])+ 
	  b[1]*(dc[0]*Pi[20] + dc[1]*Pi[21] + dc[2]*Pi[22] + dc[3]*Pi[23])+
	  b[2]*(dc[0]*Pi[24] + dc[1]*Pi[25] + dc[2]*Pi[26] + dc[3]*Pi[27])+
	  b[3]*(dc[0]*Pi[28] + dc[1]*Pi[29] + dc[2]*Pi[30] + dc[3]*Pi[31]))+
    a[2]*(b[0]*(dc[0]*Pi[32] + dc[1]*Pi[33] + dc[2]*Pi[34] + dc[3]*Pi[35])+
	  b[1]*(dc[0]*Pi[36] + dc[1]*Pi[37] + dc[2]*Pi[38] + dc[3]*Pi[39])+
	  b[2]*(dc[0]*Pi[40] + dc[1]*Pi[41] + dc[2]*Pi[42] + dc[3]*Pi[43])+
	  b[3]*(dc[0]*Pi[44] + dc[1]*Pi[45] + dc[2]*Pi[46] + dc[3]*Pi[47]))+
    a[3]*(b[0]*(dc[0]*Pi[48] + dc[1]*Pi[49] + dc[2]*Pi[50] + dc[3]*Pi[51])+
	  b[1]*(dc[0]*Pi[52] + dc[1]*Pi[53] + dc[2]*Pi[54] + dc[3]*Pi[55])+
	  b[2]*(dc[0]*Pi[56] + dc[1]*Pi[57] + dc[2]*Pi[58] + dc[3]*Pi[59])+
	  b[3]*(dc[0]*Pi[60] + dc[1]*Pi[61] + dc[2]*Pi[62] + dc[3]*Pi[63]));
}



#endif
