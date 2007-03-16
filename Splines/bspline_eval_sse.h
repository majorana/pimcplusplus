#ifndef BSPLINE_EVAL_SSE_H
#define BSPLINE_EVAL_SSE_H

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <stdio.h>
#include <math.h>

extern __m128   A0,   A1,   A2,   A3;
extern __m128  dA0,  dA1,  dA2,  dA3;
extern __m128 d2A0, d2A1, d2A2, d2A3;




/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_s (UBspline_1d_s * restrict spline, 
		    double x, float* restrict val)
{

}

/* Value and first derivative */
inline void
eval_UBspline_1d_s_vg (UBspline_1d_s * restrict spline, double x, 
		     float* restrict val, float* restrict grad)
{

}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_s_vgl (UBspline_1d_s * restrict spline, double x, 
			float* restrict val, float* restrict grad,
			float* restrict lapl)
{

}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
		    double x, double y, float* restrict val)
{

}


/* Value and gradient */
inline void
eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
		       double x, double y, 
		       float* restrict val, float* restrict grad)
{

}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_s_vgl (UBspline_2d_s * restrict spline, 
			double x, double y, float* restrict val, 
			float* restrict grad, float* restrict lapl)
{

}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_s_vgh (UBspline_2d_s * restrict spline, 
			double x, double y, float* restrict val, 
			float* restrict grad, float* restrict hess)
{

}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_s (UBspline_3d_s * restrict spline, 
		    double x, double y, double z,
		    float* restrict val)
{

}

/* Value and gradient */
inline void
eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad)
{

}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_s_vgl (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad, float* restrict lapl)
{

}




/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_s_vgh (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad, 
			float* restrict hess)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;

//    __m128 A0 = _mm_set_ps (-1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0);
//    __m128 A1 = _mm_set_ps ( 3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0);
//    __m128 A2 = _mm_set_ps (-3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0);
//    __m128 A3 = _mm_set_ps ( 1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  _mm_prefetch ((void*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,1), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,3), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,1), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,3), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,1), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,3), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,1), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,3), _MM_HINT_T0);

  __m128 tpx, tpy, tpz, a, b, c, da, db, dc, d2a, d2b, d2c,
    cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP, d2bcP, dbdcP, bd2cP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  tpx = _mm_set_ps (tx*tx*tx, tx*tx, tx, 1.0);
  tpy = _mm_set_ps (ty*ty*ty, ty*ty, ty, 1.0);
  tpz = _mm_set_ps (tz*tz*tz, tz*tz, tz, 1.0);
  
  // x-dependent vectors
  tmp0 = _mm_mul_ps (A0, tpx);
  tmp1 = _mm_mul_ps (A1, tpx);
  tmp2 = _mm_mul_ps (A2, tpx);
  tmp3 = _mm_mul_ps (A3, tpx);
  a    = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  tmp0 = _mm_mul_ps (dA0, tpx);
  tmp1 = _mm_mul_ps (dA1, tpx);
  tmp2 = _mm_mul_ps (dA2, tpx);
  tmp3 = _mm_mul_ps (dA3, tpx);
  da   = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  tmp0 = _mm_mul_ps (d2A0, tpx);
  tmp1 = _mm_mul_ps (d2A1, tpx);
  tmp2 = _mm_mul_ps (d2A2, tpx);
  tmp3 = _mm_mul_ps (d2A3, tpx);
  d2a  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  // y-dependent vectors
  tmp0 = _mm_mul_ps (A0, tpy);
  tmp1 = _mm_mul_ps (A1, tpy);
  tmp2 = _mm_mul_ps (A2, tpy);
  tmp3 = _mm_mul_ps (A3, tpy);
  b    = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  tmp0 = _mm_mul_ps (dA0, tpy);
  tmp1 = _mm_mul_ps (dA1, tpy);
  tmp2 = _mm_mul_ps (dA2, tpy);
  tmp3 = _mm_mul_ps (dA3, tpy);
  db   = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  tmp0 = _mm_mul_ps (d2A0, tpy);
  tmp1 = _mm_mul_ps (d2A1, tpy);
  tmp2 = _mm_mul_ps (d2A2, tpy);
  tmp3 = _mm_mul_ps (d2A3, tpy);
  d2b  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  // z-dependent vectors
  tmp0 = _mm_mul_ps (A0, tpz);
  tmp1 = _mm_mul_ps (A1, tpz);
  tmp2 = _mm_mul_ps (A2, tpz);
  tmp3 = _mm_mul_ps (A3, tpz);
  c    = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  tmp0 = _mm_mul_ps (dA0, tpz);
  tmp1 = _mm_mul_ps (dA1, tpz);
  tmp2 = _mm_mul_ps (dA2, tpz);
  tmp3 = _mm_mul_ps (dA3, tpz);
  dc   = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  tmp0 = _mm_mul_ps (d2A0, tpz);
  tmp1 = _mm_mul_ps (d2A1, tpz);
  tmp2 = _mm_mul_ps (d2A2, tpz);
  tmp3 = _mm_mul_ps (d2A3, tpz);
  d2c  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  // Compute cP, dcP, and d2cP products
  tmp0 = _mm_loadu_ps (P(0,0));
  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));
  tmp3 = _mm_loadu_ps (P(0,3));
  tmp4 = _mm_mul_ps (tmp0, c);
  tmp5 = _mm_mul_ps (tmp1, c);
  tmp6 = _mm_mul_ps (tmp2, c);
  tmp7 = _mm_mul_ps (tmp3, c);
  cP[0] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, dc);
  tmp5 = _mm_mul_ps (tmp1, dc);
  tmp6 = _mm_mul_ps (tmp2, dc);
  tmp7 = _mm_mul_ps (tmp3, dc);
  dcP[0] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, d2c);
  tmp5 = _mm_mul_ps (tmp1, d2c);
  tmp6 = _mm_mul_ps (tmp2, d2c);
  tmp7 = _mm_mul_ps (tmp3, d2c);
  d2cP[0] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));

  tmp0 = _mm_loadu_ps (P(1,0));
  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));
  tmp3 = _mm_loadu_ps (P(1,3));
  tmp4 = _mm_mul_ps (tmp0, c);
  tmp5 = _mm_mul_ps (tmp1, c);
  tmp6 = _mm_mul_ps (tmp2, c);
  tmp7 = _mm_mul_ps (tmp3, c);
  cP[1] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, dc);
  tmp5 = _mm_mul_ps (tmp1, dc);
  tmp6 = _mm_mul_ps (tmp2, dc);
  tmp7 = _mm_mul_ps (tmp3, dc);
  dcP[1] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, d2c);
  tmp5 = _mm_mul_ps (tmp1, d2c);
  tmp6 = _mm_mul_ps (tmp2, d2c);
  tmp7 = _mm_mul_ps (tmp3, d2c);
  d2cP[1] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));

  tmp0 = _mm_loadu_ps (P(2,0));
  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));
  tmp3 = _mm_loadu_ps (P(2,3));
  tmp4 = _mm_mul_ps (tmp0, c);
  tmp5 = _mm_mul_ps (tmp1, c);
  tmp6 = _mm_mul_ps (tmp2, c);
  tmp7 = _mm_mul_ps (tmp3, c);
  cP[2] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, dc);
  tmp5 = _mm_mul_ps (tmp1, dc);
  tmp6 = _mm_mul_ps (tmp2, dc);
  tmp7 = _mm_mul_ps (tmp3, dc);
  dcP[2] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, d2c);
  tmp5 = _mm_mul_ps (tmp1, d2c);
  tmp6 = _mm_mul_ps (tmp2, d2c);
  tmp7 = _mm_mul_ps (tmp3, d2c);
  d2cP[2] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));

  tmp0 = _mm_loadu_ps (P(3,0));
  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));
  tmp3 = _mm_loadu_ps (P(3,3));
  tmp4 = _mm_mul_ps (tmp0, c);
  tmp5 = _mm_mul_ps (tmp1, c);
  tmp6 = _mm_mul_ps (tmp2, c);
  tmp7 = _mm_mul_ps (tmp3, c);
  cP[3] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, dc);
  tmp5 = _mm_mul_ps (tmp1, dc);
  tmp6 = _mm_mul_ps (tmp2, dc);
  tmp7 = _mm_mul_ps (tmp3, dc);
  dcP[3] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  tmp4 = _mm_mul_ps (tmp0, d2c);
  tmp5 = _mm_mul_ps (tmp1, d2c);
  tmp6 = _mm_mul_ps (tmp2, d2c);
  tmp7 = _mm_mul_ps (tmp3, d2c);
  d2cP[3] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  
  tmp0 = _mm_mul_ps (b, cP[0]);
  tmp1 = _mm_mul_ps (b, cP[1]);
  tmp2 = _mm_mul_ps (b, cP[2]);
  tmp3 = _mm_mul_ps (b, cP[3]);
  bcP  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  tmp0 = _mm_mul_ps (db, cP[0]);
  tmp1 = _mm_mul_ps (db, cP[1]);
  tmp2 = _mm_mul_ps (db, cP[2]);
  tmp3 = _mm_mul_ps (db, cP[3]);
  dbcP = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  tmp0 = _mm_mul_ps (b, dcP[0]);
  tmp1 = _mm_mul_ps (b, dcP[1]);
  tmp2 = _mm_mul_ps (b, dcP[2]);
  tmp3 = _mm_mul_ps (b, dcP[3]);
  bdcP  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  tmp0 = _mm_mul_ps (d2b, cP[0]);
  tmp1 = _mm_mul_ps (d2b, cP[1]);
  tmp2 = _mm_mul_ps (d2b, cP[2]);
  tmp3 = _mm_mul_ps (d2b, cP[3]);
  d2bcP= _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  tmp0 = _mm_mul_ps (b, d2cP[0]);
  tmp1 = _mm_mul_ps (b, d2cP[1]);
  tmp2 = _mm_mul_ps (b, d2cP[2]);
  tmp3 = _mm_mul_ps (b, d2cP[3]);
  bd2cP= _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  tmp0  = _mm_mul_ps (db, dcP[0]);
  tmp1  = _mm_mul_ps (db, dcP[1]);
  tmp2  = _mm_mul_ps (db, dcP[2]);
  tmp3  = _mm_mul_ps (db, dcP[3]);
  dbdcP = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  // Compute value
  tmp0 = _mm_mul_ps (a, bcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (val, tmp0);

  // Compute gradient
  // x
  tmp0 = _mm_mul_ps (da, bcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(grad[0]), tmp0);
  // y
  tmp0 = _mm_mul_ps (a, dbcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(grad[1]), tmp0);
  // z
  tmp0 = _mm_mul_ps (a, bdcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(grad[2]), tmp0);

  // Compute hessian
  // d2x
  tmp0 = _mm_mul_ps  (d2a, bcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(hess[0]), tmp0);
  
  // d2y
  tmp0 = _mm_mul_ps  (a, d2bcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(hess[4]), tmp0);
  
  // d2z
  tmp0 = _mm_mul_ps  (a, bd2cP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(hess[8]), tmp0);

  // dx dy
  tmp0 = _mm_mul_ps  (da, dbcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(hess[1]), tmp0);

  // dx dz
  tmp0 = _mm_mul_ps  (da, bdcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(hess[2]), tmp0);

  // dy dz
  tmp0 = _mm_mul_ps  (da, dbdcP);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  tmp0 = _mm_hadd_ps (tmp0, tmp0);
  _mm_store_ss (&(hess[5]), tmp0);

  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dxInv;
  grad[2] *= dxInv;
  hess[0] *= dxInv*dxInv;
  hess[4] *= dyInv*dyInv;
  hess[8] *= dzInv*dzInv;
  hess[1] *= dxInv*dyInv;
  hess[2] *= dxInv*dzInv;
  hess[5] *= dyInv*dzInv;
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];

#undef P


}
// inline void
// find_x (float x, int* restrict index, 
// 	__m128* restrict _tp)
// {
//   __m128 _fi =_mm_set_ss (fi);
//   int i = _mm_cvttsd_si32 (_fi);
//   float t = fi - (float)i;
//   _tp = _mm_set_ps (t*t*t, t*t, t, 1.0);
// }

// inline void
// find_xy (float x, float y,
// 	 __m128 _start, __m128 _delta_inv,
// 	 int* restrict ix, int* restrict iy,
// 	 __m128* restrict tpx, __m128* restrict tpy)
// {
//   __m128 _r = _mm_set_ps (x, x, y, y);
//   __m128 _x, _y;
//   _r = _mm_sub_ps (_r, _start);
//   _r = _mm_mul_ps (_r, _delta_inv);
//   //  __m64 _ixiy = _mm_cvttps_pi32 (_r);
//   *ix = _mm_cvttss_si32 (_r);
//   __m128 _  = _mm_unpack_hi (_r, _r);
//   *iy = _mm_cvttss_si32 (_r);
// }



#endif
