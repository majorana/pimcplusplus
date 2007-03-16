#ifndef BSPLINE_EVAL_SSE_H
#define BSPLINE_EVAL_SSE_H

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>

extern __m128   _A0,   _A1,   _A2,   _A3;
extern __m128  _dA0,  _dA1,  _dA2,  _dA3;
extern __m128 _d2A0, _d2A1, _d2A2, _d2A3;

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

inline void
find_x (float x, int* restrict index, 
	__m128* restrict _tp)
{
  __m128 _fi =_mm_set_ss (fi);
  int i = _mm_cvttsd_si32 (_fi);
  float t = fi - (float)i;
  _tp = _mm_set_ps (t*t*t, t*t, t, 1.0);
}

inline void
find_xy (float x, float y,
	 __m128 _start, __m128 _delta_inv,
	 int* restrict ix, int* restrict iy,
	 __m128* restrict tpx, __m128* restrict tpy)
{
  __m128 _r = _mm_set_ps (x, x, y, y);
  __m128 _x, _y;
  _r = _mm_sub_ps (_r, _start);
  _r = _mm_mul_ps (_r, _delta_inv);
  //  __m64 _ixiy = _mm_cvttps_pi32 (_r);
  *ix = _mm_cvttss_si32 (_r);
  __m128 _  = _mm_unpack_hi (_r, _r);
  *iy = _mm_cvttss_si32 (_r);
  
  

}


/* Value only */
inline void
eval_UBspline_1d_s (UBspline_1d_s * restrict spline, 
		    double x, float* restrict val)
{
  x -= spline->x_grid.start;
  double fi = x * spline->x_grid.delta_inv;
  __m128d _fi = _mm_set_sd(fi);
  int i = _mm_cvttsd_si32 (_fi);
  double t = fi - (double)i;
  
  __m128 tp;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

#endif
