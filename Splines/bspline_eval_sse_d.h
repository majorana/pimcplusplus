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

extern __m128d 
    A0_01,   A0_23,   A1_01,   A1_23,   A2_01,   A2_23,   A3_01,   A3_23, 
   dA0_01,  dA0_23,  dA1_01,  dA1_23,  dA2_01,  dA2_23,  dA3_01,  dA3_23, 
  d2A0_01, d2A0_23, d2A1_01, d2A1_23, d2A2_01, d2A2_23, d2A3_01, d2A3_23;



/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_d (UBspline_1d_d * restrict spline, 
		    double x, double* restrict val)
{

}

/* Value and first derivative */
inline void
eval_UBspline_1d_d_vg (UBspline_1d_d * restrict spline, double x, 
		     double* restrict val, double* restrict grad)
{

}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_d_vgl (UBspline_1d_d * restrict spline, double x, 
			double* restrict val, double* restrict grad,
			double* restrict lapl)
{

}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_d (UBspline_2d_d * restrict spline, 
		    double x, double y, double* restrict val)
{

}


/* Value and gradient */
inline void
eval_UBspline_2d_d_vg (UBspline_2d_d * restrict spline, 
		       double x, double y, 
		       double* restrict val, double* restrict grad)
{

}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_d_vgl (UBspline_2d_d * restrict spline, 
			double x, double y, double* restrict val, 
			double* restrict grad, double* restrict lapl)
{

}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_d_vgh (UBspline_2d_d * restrict spline, 
			double x, double y, double* restrict val, 
			double* restrict grad, double* restrict hess)
{

}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_d (UBspline_3d_d * restrict spline, 
		    double x, double y, double z,
		    double* restrict val)
{

}

/* Value and gradient */
inline void
eval_UBspline_3d_d_vg (UBspline_3d_d * restrict spline, 
			double x, double y, double z,
			double* restrict val, double* restrict grad)
{

}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_d_vgl (UBspline_3d_d * restrict spline, 
			double x, double y, double z,
			double* restrict val, double* restrict grad, double* restrict lapl)
{

}




/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_d_vgh (UBspline_3d_d * restrict spline, 
			double x, double y, double z,
			double* restrict val, double* restrict grad, 
			double* restrict hess)
{
  _mm_prefetch ((void*)  &A0_01,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A0_23,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A1_01,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A1_23,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A2_01,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A2_23,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A3_01,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A3_23,_MM_HINT_T0);  

  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((void*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((void*)P(3,3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
    a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
    a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
    cP[8], dcP[8], d2cP[8], 
    bcP01, dbcP01, bdcP01, d2bcP01, dbdcP01, bd2cP01,
    bcP23, dbcP23, bdcP23, d2bcP23, dbdcP23, bd2cP23,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
  tpz23 = _mm_set_pd (tz, 1.0);

  
  // x-dependent vectors
  tmp0 = _mm_hadd_pd(_mm_mul_pd (A0_01, tpx01), _mm_mul_pd (A0_23, tpx23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (A1_01, tpx01), _mm_mul_pd (A1_23, tpx23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (A2_01, tpx01), _mm_mul_pd (A2_23, tpx23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (A3_01, tpx01), _mm_mul_pd (A3_23, tpx23));
  a01  = _mm_hadd_pd(tmp0, tmp1);
  a23  = _mm_hadd_pd(tmp2, tmp3);

  tmp0 = _mm_hadd_pd(_mm_mul_pd (dA0_01, tpx01), _mm_mul_pd (dA0_23, tpx23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (dA1_01, tpx01), _mm_mul_pd (dA1_23, tpx23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (dA2_01, tpx01), _mm_mul_pd (dA2_23, tpx23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (dA3_01, tpx01), _mm_mul_pd (dA3_23, tpx23));
  da01  = _mm_hadd_pd(tmp0, tmp1);
  da23  = _mm_hadd_pd(tmp2, tmp3);

  tmp0 = _mm_hadd_pd(_mm_mul_pd (d2A0_01, tpx01), _mm_mul_pd (d2A0_23, tpx23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (d2A1_01, tpx01), _mm_mul_pd (d2A1_23, tpx23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (d2A2_01, tpx01), _mm_mul_pd (d2A2_23, tpx23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (d2A3_01, tpx01), _mm_mul_pd (d2A3_23, tpx23));
  d2a01  = _mm_hadd_pd(tmp0, tmp1);
  d2a23  = _mm_hadd_pd(tmp2, tmp3);

  // y-dependent vectors
  tmp0 = _mm_hadd_pd(_mm_mul_pd (A0_01, tpy01), _mm_mul_pd (A0_23, tpy23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (A1_01, tpy01), _mm_mul_pd (A1_23, tpy23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (A2_01, tpy01), _mm_mul_pd (A2_23, tpy23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (A3_01, tpy01), _mm_mul_pd (A3_23, tpy23));
  b01  = _mm_hadd_pd(tmp0, tmp1);
  b23  = _mm_hadd_pd(tmp2, tmp3);

  tmp0 = _mm_hadd_pd(_mm_mul_pd (dA0_01, tpy01), _mm_mul_pd (dA0_23, tpy23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (dA1_01, tpy01), _mm_mul_pd (dA1_23, tpy23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (dA2_01, tpy01), _mm_mul_pd (dA2_23, tpy23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (dA3_01, tpy01), _mm_mul_pd (dA3_23, tpy23));
  db01  = _mm_hadd_pd(tmp0, tmp1);
  db23  = _mm_hadd_pd(tmp2, tmp3);

  tmp0 = _mm_hadd_pd(_mm_mul_pd (d2A0_01, tpy01), _mm_mul_pd (d2A0_23, tpz23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (d2A1_01, tpy01), _mm_mul_pd (d2A1_23, tpz23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (d2A2_01, tpy01), _mm_mul_pd (d2A2_23, tpz23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (d2A3_01, tpy01), _mm_mul_pd (d2A3_23, tpz23));
  d2b01  = _mm_hadd_pd(tmp0, tmp1);
  d2b23  = _mm_hadd_pd(tmp2, tmp3);

  // z-dependent vectors
  tmp0 = _mm_hadd_pd(_mm_mul_pd (A0_01, tpz01), _mm_mul_pd (A0_23, tpz23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (A1_01, tpz01), _mm_mul_pd (A1_23, tpz23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (A2_01, tpz01), _mm_mul_pd (A2_23, tpz23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (A3_01, tpz01), _mm_mul_pd (A3_23, tpz23));
  c01  = _mm_hadd_pd(tmp0, tmp1);
  c23  = _mm_hadd_pd(tmp2, tmp3);

  tmp0 = _mm_hadd_pd(_mm_mul_pd (dA0_01, tpz01), _mm_mul_pd (dA0_23, tpz23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (dA1_01, tpz01), _mm_mul_pd (dA1_23, tpz23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (dA2_01, tpz01), _mm_mul_pd (dA2_23, tpz23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (dA3_01, tpz01), _mm_mul_pd (dA3_23, tpz23));
  dc01  = _mm_hadd_pd(tmp0, tmp1);
  dc23  = _mm_hadd_pd(tmp2, tmp3);

  tmp0 = _mm_hadd_pd(_mm_mul_pd (d2A0_01, tpz01), _mm_mul_pd (d2A0_23, tpz23));
  tmp1 = _mm_hadd_pd(_mm_mul_pd (d2A1_01, tpz01), _mm_mul_pd (d2A1_23, tpz23));
  tmp2 = _mm_hadd_pd(_mm_mul_pd (d2A2_01, tpz01), _mm_mul_pd (d2A2_23, tpz23));
  tmp3 = _mm_hadd_pd(_mm_mul_pd (d2A3_01, tpz01), _mm_mul_pd (d2A3_23, tpz23));
  d2c01  = _mm_hadd_pd(tmp0, tmp1);
  d2c23  = _mm_hadd_pd(tmp2, tmp3);
 
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st eighth
  tmp0    = _mm_loadu_pd (P(0,0,0));
  tmp1    = _mm_loadu_pd (P(0,0,2));
  tmp2    = _mm_loadu_pd (P(0,1,0));
  tmp3    = _mm_loadu_pd (P(0,1,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[0]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[0]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[0] = _mm_hadd_pd (tmp4, tmp5);

  // 2nd eighth
  tmp0    = _mm_loadu_pd (P(0,2,0));
  tmp1    = _mm_loadu_pd (P(0,2,2));
  tmp2    = _mm_loadu_pd (P(0,3,0));
  tmp3    = _mm_loadu_pd (P(0,3,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[1]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[1]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[1] = _mm_hadd_pd (tmp4, tmp5);

  // 3rd eighth
  tmp0    = _mm_loadu_pd (P(1,0,0));
  tmp1    = _mm_loadu_pd (P(1,0,2));
  tmp2    = _mm_loadu_pd (P(1,1,0));
  tmp3    = _mm_loadu_pd (P(1,1,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[2]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[2]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[2] = _mm_hadd_pd (tmp4, tmp5);

  // 4th eighth
  tmp0    = _mm_loadu_pd (P(1,2,0));
  tmp1    = _mm_loadu_pd (P(1,2,2));
  tmp2    = _mm_loadu_pd (P(1,3,0));
  tmp3    = _mm_loadu_pd (P(1,3,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[3]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[3]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[3] = _mm_hadd_pd (tmp4, tmp5);

  // 5th eighth
  tmp0    = _mm_loadu_pd (P(2,0,0));
  tmp1    = _mm_loadu_pd (P(2,0,2));
  tmp2    = _mm_loadu_pd (P(2,1,0));
  tmp3    = _mm_loadu_pd (P(2,1,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[4]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[4]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[4] = _mm_hadd_pd (tmp4, tmp5);

  // 6th eighth
  tmp0    = _mm_loadu_pd (P(2,2,0));
  tmp1    = _mm_loadu_pd (P(2,2,2));
  tmp2    = _mm_loadu_pd (P(2,3,0));
  tmp3    = _mm_loadu_pd (P(2,3,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[5]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[5]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[5] = _mm_hadd_pd (tmp4, tmp5);

  // 7th eighth
  tmp0    = _mm_loadu_pd (P(3,0,0));
  tmp1    = _mm_loadu_pd (P(3,0,2));
  tmp2    = _mm_loadu_pd (P(3,1,0));
  tmp3    = _mm_loadu_pd (P(3,1,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[6]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[6]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[6] = _mm_hadd_pd (tmp4, tmp5);

  // 8th eighth
  tmp0    = _mm_loadu_pd (P(3,2,0));
  tmp1    = _mm_loadu_pd (P(3,2,2));
  tmp2    = _mm_loadu_pd (P(3,3,0));
  tmp3    = _mm_loadu_pd (P(3,3,2));
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, c01),   _mm_mul_pd (tmp1, c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, c01),   _mm_mul_pd (tmp3, c23));
  cP[7]   = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, dc01),  _mm_mul_pd (tmp1, dc23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, dc01),  _mm_mul_pd (tmp3, dc23));
  dcP[7]  = _mm_hadd_pd (tmp4, tmp5);
  tmp4    = _mm_hadd_pd (_mm_mul_pd (tmp0, d2c01), _mm_mul_pd (tmp1, d2c23));
  tmp5    = _mm_hadd_pd (_mm_mul_pd (tmp2, d2c01), _mm_mul_pd (tmp3, d2c23));
  d2cP[7] = _mm_hadd_pd (tmp4, tmp5);

  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  tmp0 = _mm_hadd_pd (_mm_mul_pd (b01, cP[0]), _mm_mul_pd (b23, cP[1]));
  tmp1 = _mm_hadd_pd (_mm_mul_pd (b01, cP[2]), _mm_mul_pd (b23, cP[3]));
  tmp2 = _mm_hadd_pd (_mm_mul_pd (b01, cP[4]), _mm_mul_pd (b23, cP[5]));
  tmp3 = _mm_hadd_pd (_mm_mul_pd (b01, cP[6]), _mm_mul_pd (b23, cP[7]));
  bcP01 = _mm_hadd_pd (tmp0, tmp1);
  bcP23 = _mm_hadd_pd (tmp2, tmp3);

  tmp0 = _mm_hadd_pd (_mm_mul_pd (db01, cP[0]), _mm_mul_pd (db23, cP[1]));
  tmp1 = _mm_hadd_pd (_mm_mul_pd (db01, cP[2]), _mm_mul_pd (db23, cP[3]));
  tmp2 = _mm_hadd_pd (_mm_mul_pd (db01, cP[4]), _mm_mul_pd (db23, cP[5]));
  tmp3 = _mm_hadd_pd (_mm_mul_pd (db01, cP[6]), _mm_mul_pd (db23, cP[7]));
  dbcP01 = _mm_hadd_pd (tmp0, tmp1);
  dbcP23 = _mm_hadd_pd (tmp2, tmp3);

  tmp0 = _mm_hadd_pd (_mm_mul_pd (b01, dcP[0]), _mm_mul_pd (b23, dcP[1]));
  tmp1 = _mm_hadd_pd (_mm_mul_pd (b01, dcP[2]), _mm_mul_pd (b23, dcP[3]));
  tmp2 = _mm_hadd_pd (_mm_mul_pd (b01, dcP[4]), _mm_mul_pd (b23, dcP[5]));
  tmp3 = _mm_hadd_pd (_mm_mul_pd (b01, dcP[6]), _mm_mul_pd (b23, dcP[7]));
  bdcP01 = _mm_hadd_pd (tmp0, tmp1);
  bdcP23 = _mm_hadd_pd (tmp2, tmp3);

  tmp0 = _mm_hadd_pd (_mm_mul_pd (d2b01, cP[0]), _mm_mul_pd (d2b23, cP[1]));
  tmp1 = _mm_hadd_pd (_mm_mul_pd (d2b01, cP[2]), _mm_mul_pd (d2b23, cP[3]));
  tmp2 = _mm_hadd_pd (_mm_mul_pd (d2b01, cP[4]), _mm_mul_pd (d2b23, cP[5]));
  tmp3 = _mm_hadd_pd (_mm_mul_pd (d2b01, cP[6]), _mm_mul_pd (d2b23, cP[7]));
  d2bcP01 = _mm_hadd_pd (tmp0, tmp1);
  d2bcP23 = _mm_hadd_pd (tmp2, tmp3);

  tmp0 = _mm_hadd_pd (_mm_mul_pd (b01, d2cP[0]), _mm_mul_pd (b23, d2cP[1]));
  tmp1 = _mm_hadd_pd (_mm_mul_pd (b01, d2cP[2]), _mm_mul_pd (b23, d2cP[3]));
  tmp2 = _mm_hadd_pd (_mm_mul_pd (b01, d2cP[4]), _mm_mul_pd (b23, d2cP[5]));
  tmp3 = _mm_hadd_pd (_mm_mul_pd (b01, d2cP[6]), _mm_mul_pd (b23, d2cP[7]));
  bd2cP01 = _mm_hadd_pd (tmp0, tmp1);
  bd2cP23 = _mm_hadd_pd (tmp2, tmp3);
  
  tmp0 = _mm_hadd_pd (_mm_mul_pd (db01, dcP[0]), _mm_mul_pd (db23, dcP[1]));
  tmp1 = _mm_hadd_pd (_mm_mul_pd (db01, dcP[2]), _mm_mul_pd (db23, dcP[3]));
  tmp2 = _mm_hadd_pd (_mm_mul_pd (db01, dcP[4]), _mm_mul_pd (db23, dcP[5]));
  tmp3 = _mm_hadd_pd (_mm_mul_pd (db01, dcP[6]), _mm_mul_pd (db23, dcP[7]));
  dbdcP01 = _mm_hadd_pd (tmp0, tmp1);
  dbdcP23 = _mm_hadd_pd (tmp2, tmp3);

  // Compute value
  tmp0 = _mm_hadd_pd (_mm_mul_pd (a01, bcP01), _mm_mul_pd (a23, bcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (val, tmp0);
  // Compute gradient
  // x
  tmp0 = _mm_hadd_pd (_mm_mul_pd (da01, bcP01), _mm_mul_pd (da23, bcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&(grad[0]), tmp0);
  // y
  tmp0 = _mm_hadd_pd (_mm_mul_pd (a01, dbcP01), _mm_mul_pd (a23, dbcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&grad[1], tmp0);
  // z
  tmp0 = _mm_hadd_pd (_mm_mul_pd (a01, bdcP01), _mm_mul_pd (a23, bdcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&grad[2], tmp0);

  // Compute hessian
  // d2x
  tmp0 = _mm_hadd_pd (_mm_mul_pd (d2a01, bcP01), _mm_mul_pd (d2a23, bcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&hess[0], tmp0);
  // d2y
  tmp0 = _mm_hadd_pd (_mm_mul_pd (a01, d2bcP01), _mm_mul_pd (a23, d2bcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&hess[4], tmp0);
  // d2z
  tmp0 = _mm_hadd_pd (_mm_mul_pd (a01, bd2cP01), _mm_mul_pd (a23, bd2cP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&hess[8], tmp0);
  // dx dy
  tmp0 = _mm_hadd_pd (_mm_mul_pd (da01, dbcP01), _mm_mul_pd (da23, dbcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&hess[1], tmp0);
  // dx dz
  tmp0 = _mm_hadd_pd (_mm_mul_pd (da01, bdcP01), _mm_mul_pd (da23, bdcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&hess[2], tmp0);
  // dy dz
  tmp0 = _mm_hadd_pd (_mm_mul_pd (a01, dbdcP01), _mm_mul_pd (a23, dbdcP23));
  tmp0 = _mm_hadd_pd (tmp0, tmp0);
  _mm_store_sd (&hess[5], tmp0);
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
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


#endif
