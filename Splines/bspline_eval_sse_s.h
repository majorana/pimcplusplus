#ifndef BSPLINE_EVAL_SSE_S_H
#define BSPLINE_EVAL_SSE_S_H

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
/// SSE3 add "horizontal add" instructions, which makes things
/// simpler and faster
#ifdef __SSE3__
#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  __m128 r0 = _mm_hadd_ps (_mm_mul_ps (M0, v), _mm_mul_ps (M1, v)); \
  __m128 r1 = _mm_hadd_ps (_mm_mul_ps (M2, v), _mm_mul_ps (M3, v)); \
  r = _mm_hadd_ps (r0, r1);                                         \
 } while (0);
#define _MM_DOT4_PS(A, B, p)                                       \
do {                                                                \
  __m128 t  = _mm_mul_ps (A, B);                                    \
  __m128 t1 = _mm_hadd_ps (t,t);                                    \
  __m128 r  = _mm_hadd_ps (t1, t1);                                 \
  _mm_store_ss (&(p), r);                                           \
} while(0);
#else
// Use plain-old SSE instructions
#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  __m128 r0 = _mm_mul_ps (M0, v);                                   \
  __m128 r1 = _mm_mul_ps (M1, v);				    \
  __m128 r2 = _mm_mul_ps (M2, v);                                   \
  __m128 r3 = _mm_mul_ps (M3, v);				    \
  _MM_TRANSPOSE4_PS (r0, r1, r2, r3);                               \
  r = _mm_add_ps (_mm_add_ps (r0, r1), _mm_add_ps (r2, r3));        \
 } while (0);
#define _MM_DOT4_PS(A, B, p)                                        \
do {                                                                \
  __m128 t    = _mm_mul_ps (A, B);                                  \
  __m128 alo  = _mm_shuffle_ps (t, t, _MM_SHUFFLE(0,1,0,1));	    \
  __m128 ahi  = _mm_shuffle_ps (t, t, _MM_SHUFFLE(2,3,2,3));	    \
  __m128 a    = _mm_add_ps (alo, ahi);                              \
  __m128 rlo  = _mm_shuffle_ps (a, a, _MM_SHUFFLE(0,0,0,0));	     \
  __m128 rhi  = _mm_shuffle_ps (a, a, _MM_SHUFFLE(1,1,1,1));	     \
  __m128 r    = _mm_add_ps (rlo, rhi);                              \
  _mm_store_ss (&(p), r);                                           \
} while(0);
#endif

/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_s_vgh (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad, 
			float* restrict hess)
{
  _mm_prefetch ((void*)  &A0,_MM_HINT_T0);  _mm_prefetch ((void*)  &A1,_MM_HINT_T0);  
  _mm_prefetch ((void*)  &A2,_MM_HINT_T0);  _mm_prefetch ((void*)  &A3,_MM_HINT_T0);

  /// SSE mesh point determination
  __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
  __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				 spline->z_grid.start, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				 spline->z_grid.delta_inv, 0.0);
  xyz = _mm_sub_ps (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
  // intpart = trunc (ux, uy, uz)
  __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
  __m128i ixiyiz;
  _mm_storeu_si128 (&ixiyiz, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
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

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txtytz, txtytz);
  __m128 t3     = _mm_mul_ps (t2, txtytz);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txtytz;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP, d2bcP, dbdcP, bd2cP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  
  // x-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpx,   a);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpx,  da);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpx, d2a);

//   tmp0 = _mm_mul_ps (A0, tpx);
//   tmp1 = _mm_mul_ps (A1, tpx);
//   tmp2 = _mm_mul_ps (A2, tpx);
//   tmp3 = _mm_mul_ps (A3, tpx);
//   a    = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

//   tmp0 = _mm_mul_ps (dA0, tpx);
//   tmp1 = _mm_mul_ps (dA1, tpx);
//   tmp2 = _mm_mul_ps (dA2, tpx);
//   tmp3 = _mm_mul_ps (dA3, tpx);
//   da   = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

//   tmp0 = _mm_mul_ps (d2A0, tpx);
//   tmp1 = _mm_mul_ps (d2A1, tpx);
//   tmp2 = _mm_mul_ps (d2A2, tpx);
//   tmp3 = _mm_mul_ps (d2A3, tpx);
//   d2a  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  // y-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpy,   b);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpy,  db);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpy, d2b);

//   tmp0 = _mm_mul_ps (A0, tpy);
//   tmp1 = _mm_mul_ps (A1, tpy);
//   tmp2 = _mm_mul_ps (A2, tpy);
//   tmp3 = _mm_mul_ps (A3, tpy);
//   b    = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

//   tmp0 = _mm_mul_ps (dA0, tpy);
//   tmp1 = _mm_mul_ps (dA1, tpy);
//   tmp2 = _mm_mul_ps (dA2, tpy);
//   tmp3 = _mm_mul_ps (dA3, tpy);
//   db   = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

//   tmp0 = _mm_mul_ps (d2A0, tpy);
//   tmp1 = _mm_mul_ps (d2A1, tpy);
//   tmp2 = _mm_mul_ps (d2A2, tpy);
//   tmp3 = _mm_mul_ps (d2A3, tpy);
//   d2b  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  // z-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpz,   c);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpz,  dc);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpz, d2c);

//   tmp0 = _mm_mul_ps (A0, tpz);
//   tmp1 = _mm_mul_ps (A1, tpz);
//   tmp2 = _mm_mul_ps (A2, tpz);
//   tmp3 = _mm_mul_ps (A3, tpz);
//   c    = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

//   tmp0 = _mm_mul_ps (dA0, tpz);
//   tmp1 = _mm_mul_ps (dA1, tpz);
//   tmp2 = _mm_mul_ps (dA2, tpz);
//   tmp3 = _mm_mul_ps (dA3, tpz);
//   dc   = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

//   tmp0 = _mm_mul_ps (d2A0, tpz);
//   tmp1 = _mm_mul_ps (d2A1, tpz);
//   tmp2 = _mm_mul_ps (d2A2, tpz);
//   tmp3 = _mm_mul_ps (d2A3, tpz);
//   d2c  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps (tmp2, tmp3));

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0));
  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));
  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[0]);
//   tmp4 = _mm_mul_ps (tmp0, c);
//   tmp5 = _mm_mul_ps (tmp1, c);
//   tmp6 = _mm_mul_ps (tmp2, c);
//   tmp7 = _mm_mul_ps (tmp3, c);
//   cP[0] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, dc);
//   tmp5 = _mm_mul_ps (tmp1, dc);
//   tmp6 = _mm_mul_ps (tmp2, dc);
//   tmp7 = _mm_mul_ps (tmp3, dc);
//   dcP[0] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, d2c);
//   tmp5 = _mm_mul_ps (tmp1, d2c);
//   tmp6 = _mm_mul_ps (tmp2, d2c);
//   tmp7 = _mm_mul_ps (tmp3, d2c);
//   d2cP[0] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));
  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));
  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[1]);
//   tmp4 = _mm_mul_ps (tmp0, c);
//   tmp5 = _mm_mul_ps (tmp1, c);
//   tmp6 = _mm_mul_ps (tmp2, c);
//   tmp7 = _mm_mul_ps (tmp3, c);
//   cP[1] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, dc);
//   tmp5 = _mm_mul_ps (tmp1, dc);
//   tmp6 = _mm_mul_ps (tmp2, dc);
//   tmp7 = _mm_mul_ps (tmp3, dc);
//   dcP[1] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, d2c);
//   tmp5 = _mm_mul_ps (tmp1, d2c);
//   tmp6 = _mm_mul_ps (tmp2, d2c);
//   tmp7 = _mm_mul_ps (tmp3, d2c);
//   d2cP[1] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));
  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));
  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[2]);
//   tmp4 = _mm_mul_ps (tmp0, c);
//   tmp5 = _mm_mul_ps (tmp1, c);
//   tmp6 = _mm_mul_ps (tmp2, c);
//   tmp7 = _mm_mul_ps (tmp3, c);
//   cP[2] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, dc);
//   tmp5 = _mm_mul_ps (tmp1, dc);
//   tmp6 = _mm_mul_ps (tmp2, dc);
//   tmp7 = _mm_mul_ps (tmp3, dc);
//   dcP[2] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, d2c);
//   tmp5 = _mm_mul_ps (tmp1, d2c);
//   tmp6 = _mm_mul_ps (tmp2, d2c);
//   tmp7 = _mm_mul_ps (tmp3, d2c);
//   d2cP[2] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));
  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));
  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[3]);
//   tmp4 = _mm_mul_ps (tmp0, c);
//   tmp5 = _mm_mul_ps (tmp1, c);
//   tmp6 = _mm_mul_ps (tmp2, c);
//   tmp7 = _mm_mul_ps (tmp3, c);
//   cP[3] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, dc);
//   tmp5 = _mm_mul_ps (tmp1, dc);
//   tmp6 = _mm_mul_ps (tmp2, dc);
//   tmp7 = _mm_mul_ps (tmp3, dc);
//   dcP[3] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
//   tmp4 = _mm_mul_ps (tmp0, d2c);
//   tmp5 = _mm_mul_ps (tmp1, d2c);
//   tmp6 = _mm_mul_ps (tmp2, d2c);
//   tmp7 = _mm_mul_ps (tmp3, d2c);
//   d2cP[3] = _mm_hadd_ps (_mm_hadd_ps (tmp4, tmp5), _mm_hadd_ps(tmp6, tmp7));
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3], d2b, d2bcP);
  _MM_MATVEC4_PS (d2cP[0], d2cP[1], d2cP[2], d2cP[3],   b, bd2cP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],  db, dbdcP);
//   tmp0 = _mm_mul_ps (b, cP[0]);
//   tmp1 = _mm_mul_ps (b, cP[1]);
//   tmp2 = _mm_mul_ps (b, cP[2]);
//   tmp3 = _mm_mul_ps (b, cP[3]);
//   bcP  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

//   tmp0 = _mm_mul_ps (db, cP[0]);
//   tmp1 = _mm_mul_ps (db, cP[1]);
//   tmp2 = _mm_mul_ps (db, cP[2]);
//   tmp3 = _mm_mul_ps (db, cP[3]);
//   dbcP = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

//   tmp0 = _mm_mul_ps (b, dcP[0]);
//   tmp1 = _mm_mul_ps (b, dcP[1]);
//   tmp2 = _mm_mul_ps (b, dcP[2]);
//   tmp3 = _mm_mul_ps (b, dcP[3]);
//   bdcP  = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

//   tmp0 = _mm_mul_ps (d2b, cP[0]);
//   tmp1 = _mm_mul_ps (d2b, cP[1]);
//   tmp2 = _mm_mul_ps (d2b, cP[2]);
//   tmp3 = _mm_mul_ps (d2b, cP[3]);
//   d2bcP= _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

//   tmp0 = _mm_mul_ps (b, d2cP[0]);
//   tmp1 = _mm_mul_ps (b, d2cP[1]);
//   tmp2 = _mm_mul_ps (b, d2cP[2]);
//   tmp3 = _mm_mul_ps (b, d2cP[3]);
//   bd2cP= _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

//   tmp0  = _mm_mul_ps (db, dcP[0]);
//   tmp1  = _mm_mul_ps (db, dcP[1]);
//   tmp2  = _mm_mul_ps (db, dcP[2]);
//   tmp3  = _mm_mul_ps (db, dcP[3]);
//   dbdcP = _mm_hadd_ps (_mm_hadd_ps (tmp0, tmp1), _mm_hadd_ps(tmp2, tmp3));

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bcP, grad[0]);
  _MM_DOT4_PS (a, dbcP, grad[1]);
  _MM_DOT4_PS (a, bdcP, grad[2]);
  // Compute hessian
  _MM_DOT4_PS (d2a, bcP, hess[0]);
  _MM_DOT4_PS (a, d2bcP, hess[4]);
  _MM_DOT4_PS (a, bd2cP, hess[8]);
  _MM_DOT4_PS (da, dbcP, hess[1]);
  _MM_DOT4_PS (da, bdcP, hess[2]);
  _MM_DOT4_PS (a, dbdcP, hess[5]);

//   tmp0 = _mm_mul_ps (a, bcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (val, tmp0);

  // x
//   tmp0 = _mm_mul_ps (da, bcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(grad[0]), tmp0);
//   // y
//   tmp0 = _mm_mul_ps (a, dbcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(grad[1]), tmp0);
//   // z
//   tmp0 = _mm_mul_ps (a, bdcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(grad[2]), tmp0);

  // d2x
//   tmp0 = _mm_mul_ps  (d2a, bcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(hess[0]), tmp0);
  
  // d2y
//   tmp0 = _mm_mul_ps  (a, d2bcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(hess[4]), tmp0);
  
  // d2z
//   tmp0 = _mm_mul_ps  (a, bd2cP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(hess[8]), tmp0);

  // dx dy
//   tmp0 = _mm_mul_ps  (da, dbcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(hess[1]), tmp0);

  // dx dz
//   tmp0 = _mm_mul_ps  (da, bdcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(hess[2]), tmp0);

  // dy dz
//   tmp0 = _mm_mul_ps  (da, dbdcP);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   tmp0 = _mm_hadd_ps (tmp0, tmp0);
//   _mm_store_ss (&(hess[5]), tmp0);

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
