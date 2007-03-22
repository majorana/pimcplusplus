#ifndef BSPLINE_EVAL_SSE_C_H
#define BSPLINE_EVAL_SSE_C_H

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <stdio.h>
#include <math.h>

extern __m128   A0,   A1,   A2,   A3;
extern __m128  dA0,  dA1,  dA2,  dA3;
extern __m128 d2A0, d2A1, d2A2, d2A3;


/// SSE3 adds "horizontal add" instructions, which makes things
/// simpler and faster
#ifdef __SSE3__
#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  __m128 g0 = _mm_hadd_ps (_mm_mul_ps (M0, v), _mm_mul_ps (M1, v)); \
  __m128 g1 = _mm_hadd_ps (_mm_mul_ps (M2, v), _mm_mul_ps (M3, v)); \
  r = _mm_hadd_ps (g0, g1);                                         \
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
  __m128 g0 = _mm_mul_ps (M0, v);                                   \
  __m128 g1 = _mm_mul_ps (M1, v);				    \
  __m128 g2 = _mm_mul_ps (M2, v);                                   \
  __m128 g3 = _mm_mul_ps (M3, v);				    \
  _MM_TRANSPOSE4_PS (g0, g1, g2, g3);                               \
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


/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_c (UBspline_1d_c * restrict spline, 
		    double x, complex_float* restrict val)
{

}

/* Value and first derivative */
inline void
eval_UBspline_1d_c_vg (UBspline_1d_c * restrict spline, double x, 
		     complex_float* restrict val, complex_float* restrict grad)
{

}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_c_vgl (UBspline_1d_c * restrict spline, double x, 
			complex_float* restrict val, 
			complex_float* restrict grad,
			complex_float* restrict lapl)
{

}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_c (UBspline_2d_c * restrict spline, 
		    double x, double y, complex_float* restrict val)
{

}


/* Value and gradient */
inline void
eval_UBspline_2d_c_vg (UBspline_2d_c * restrict spline, 
		       double x, double y, 
		       complex_float* restrict val, complex_float* restrict grad)
{

}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_c_vgl (UBspline_2d_c * restrict spline, 
			double x, double y, complex_float* restrict val, 
			complex_float* restrict grad, complex_float* restrict lapl)
{

}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_c_vgh (UBspline_2d_c * restrict spline, 
			double x, double y, complex_float* restrict val, 
			complex_float* restrict grad, 
			complex_float* restrict hess)
{

}




/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_c (UBspline_3d_c * restrict spline, 
		    double x, double y, double z,
		    complex_float* restrict val)
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
  __m128 a, b, c, cP[4],bcP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  // x-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpx,   a);
  // y-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpy,   b);
  // z-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpz,   c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);

#undef P
}

/* Value and gradient */
inline void
eval_UBspline_3d_c_vg (UBspline_3d_c * restrict spline, 
		       double x, double y, double z,
		       complex_float* restrict val, 
		       complex_float* restrict grad)
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
  __m128 a, b, c, da, db, dc,
    cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  
  // x-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpx,   a);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpx,  da);
  // y-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpy,   b);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpy,  db);
  // z-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpz,   c);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpz,  dc);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bcP, grad[0]);
  _MM_DOT4_PS (a, dbcP, grad[1]);
  _MM_DOT4_PS (a, bdcP, grad[2]);
  
  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
#undef P
}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_c_vgl (UBspline_3d_c * restrict spline, 
			double x, double y, double z,
			complex_float* restrict val, 
			complex_float* restrict grad, 
			complex_float* restrict lapl)
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
  // y-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpy,   b);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpy,  db);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpy, d2b);
  // z-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpz,   c);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpz,  dc);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpz, d2c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));
  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));
  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3], d2b, d2bcP);
  _MM_MATVEC4_PS (d2cP[0], d2cP[1], d2cP[2], d2cP[3],   b, bd2cP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],  db, dbdcP);

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bcP, grad[0]);
  _MM_DOT4_PS (a, dbcP, grad[1]);
  _MM_DOT4_PS (a, bdcP, grad[2]);
  // Compute laplacian
  float lx, ly, lz;
  _MM_DOT4_PS (d2a, bcP, lx);
  _MM_DOT4_PS (a, d2bcP, ly);
  _MM_DOT4_PS (a, bd2cP, lz);
  
  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  lx *= dxInv*dxInv;
  ly *= dyInv*dyInv;
  lz *= dzInv*dzInv;
  *lapl = lx + ly + lz;	       
#undef P
}


inline void
print__m128 (__m128 val)
{
  float v[4];
  __m128 vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(0,0,0,0));
  _mm_store_ss (&(v[0]), vshuf);
  vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(1,1,1,1));
  _mm_store_ss (&(v[1]), vshuf);
  vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(2,2,2,2));
  _mm_store_ss (&(v[2]), vshuf);
  vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(3,3,3,3));
  _mm_store_ss (&(v[3]), vshuf);
  
  fprintf (stderr, "[ %8.5f, %8.5f, %8.5f, %8.5f ]", v[0], v[1], v[2], v[3]);
}


/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_c_vgh (UBspline_3d_c * restrict spline, 
			double x, double y, double z,
			complex_float* restrict val, 
			complex_float* restrict grad, 
			complex_float* restrict hess)
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
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
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
    cPr[4], dcPr[4], d2cPr[4], bcPr, dbcPr, bdcPr, d2bcPr, dbdcPr, bd2cPr,
    cPi[4], dcPi[4], d2cPi[4], bcPi, dbcPi, bdcPi, d2bcPi, dbdcPi, bd2cPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  // x-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpx,   a);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpx,  da);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpy,   b);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpy,  db);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpy, d2b);
  // z-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpz,   c);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpz,  dc);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpz, d2c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0,0));  
  tmp1 = _mm_loadu_ps (P(0,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
//   fprintf (stderr, "P(0,0,0) = [ %8.5f, %8.5f, %8.5f, %8.5f]\n",
// 	   *P(0,0,0), *(P(0,0,0)+1), *(P(0,0,0)+2), *(P(0,0,0)+3));
//   fprintf (stderr, "P(0,0,2) = [ %8.5f, %8.5f, %8.5f, %8.5f]\n",
// 	   *P(0,0,2), *(P(0,0,2)+1), *(P(0,0,2)+2), *(P(0,0,2)+3));

//   fprintf (stderr, "tmp0 = "); print__m128 (tmp0); fprintf (stderr, "\n");
//   fprintf (stderr, "tmp1 = "); print__m128 (tmp1); fprintf (stderr, "\n");
//   fprintf (stderr, "r0   = "); print__m128 (r0);   fprintf (stderr, "\n");
//   fprintf (stderr, "i0   = "); print__m128 (i0);   fprintf (stderr, "\n");
  tmp0 = _mm_loadu_ps (P(0,1,0));  tmp1 = _mm_loadu_ps (P(0,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,2,0));  tmp1 = _mm_loadu_ps (P(0,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,3,0));  tmp1 = _mm_loadu_ps (P(0,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0,0));  tmp1 = _mm_loadu_ps (P(1,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,1,0));  tmp1 = _mm_loadu_ps (P(1,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,2,0));  tmp1 = _mm_loadu_ps (P(1,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,3,0));  tmp1 = _mm_loadu_ps (P(1,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0,0));  tmp1 = _mm_loadu_ps (P(2,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,1,0));  tmp1 = _mm_loadu_ps (P(2,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,2,0));  tmp1 = _mm_loadu_ps (P(2,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,3,0));  tmp1 = _mm_loadu_ps (P(2,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0,0));  tmp1 = _mm_loadu_ps (P(3,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,1,0));  tmp1 = _mm_loadu_ps (P(3,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,2,0));  tmp1 = _mm_loadu_ps (P(3,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,3,0));  tmp1 = _mm_loadu_ps (P(3,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],   b,   bcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],  db,  dbcPr);
  _MM_MATVEC4_PS ( dcPr[0],  dcPr[1],  dcPr[2],  dcPr[3],   b,  bdcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3], d2b, d2bcPr);
  _MM_MATVEC4_PS (d2cPr[0], d2cPr[1], d2cPr[2], d2cPr[3],   b, bd2cPr);
  _MM_MATVEC4_PS ( dcPr[0],  dcPr[1],  dcPr[2],  dcPr[3],  db, dbdcPr);

  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],   b,   bcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],  db,  dbcPi);
  _MM_MATVEC4_PS ( dcPi[0],  dcPi[1],  dcPi[2],  dcPi[3],   b,  bdcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3], d2b, d2bcPi);
  _MM_MATVEC4_PS (d2cPi[0], d2cPi[1], d2cPi[2], d2cPi[3],   b, bd2cPi);
  _MM_MATVEC4_PS ( dcPi[0],  dcPi[1],  dcPi[2],  dcPi[3],  db, dbdcPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;
  float *gradr2 = ((float *)grad)+4;
  float *gradi2 = ((float *)grad)+5;
  
  // Compute value
  _MM_DOT4_PS (a, bcPr, *valr);
  _MM_DOT4_PS (a, bcPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bcPr, *gradr0);
  _MM_DOT4_PS (a, dbcPr, *gradr1);
  _MM_DOT4_PS (a, bdcPr, *gradr2);
  _MM_DOT4_PS (da, bcPi, *gradi0);
  _MM_DOT4_PS (a, dbcPi, *gradi1);
  _MM_DOT4_PS (a, bdcPi, *gradi2);
  // Compute hessian
  _MM_DOT4_PS (d2a, bcPr, *(float*)(&hess[0]));
  _MM_DOT4_PS (a, d2bcPr, *(float*)(&hess[4]));
  _MM_DOT4_PS (a, bd2cPr, *(float*)(&hess[8]));
  _MM_DOT4_PS (da, dbcPr, *(float*)(&hess[1]));
  _MM_DOT4_PS (da, bdcPr, *(float*)(&hess[2]));
  _MM_DOT4_PS (a, dbdcPr, *(float*)(&hess[5]));

  _MM_DOT4_PS (d2a, bcPi, *((float*)(&hess[0])+1));
  _MM_DOT4_PS (a, d2bcPi, *((float*)(&hess[4])+1));
  _MM_DOT4_PS (a, bd2cPi, *((float*)(&hess[8])+1));
  _MM_DOT4_PS (da, dbcPi, *((float*)(&hess[1])+1));
  _MM_DOT4_PS (da, bdcPi, *((float*)(&hess[2])+1));
  _MM_DOT4_PS (a, dbdcPi, *((float*)(&hess[5])+1));


  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
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
