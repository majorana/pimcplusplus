#ifndef BSPLINE_H
#define BSPLINE_H

// Conventions:
// Postfixes:  
// s:  single precision real
// d:  double precision real
// c:  single precision complex
// z:  double precision complex

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Basic type declarations               ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

typedef enum { PERIODIC, DERIV1, DERIV2, FLAT, NATURAL } bc_code;
typedef enum { U1D, U2D, U3D, NU1D, NU2D, NU3D } spline_code;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal, rVal;
} BCtype_s;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal, rVal;
} BCtype_d;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal_r, lVal_i, rVal_r, rVal_i;
} BCtype_c;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal_r, lVal_i, rVal_r, rVal_i;
} BCtype_z;


typedef struct
{
  double start, end;
  int num;

  // private
  double delta, delta_inv;
} Ugrid;

typedef struct
{
  spline_code code;
} Bspline;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////           Bspline structure definitions            ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
#ifdef __SSE3__
#include "bspline_structs_sse.h"
#include "bspline_eval_sse_s.h"
#include "bspline_eval_sse_d.h"
#else
#include "bspline_structs_std.h"
#include "bspline_eval_std_s.h"
#include "bspline_eval_std_d.h"
#endif


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Spline creation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/////////////////////////////////////
// Uniform, single precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
UBspline_1d_s *
create_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, float *data);

// Create 2D uniform single-precision, real Bspline
UBspline_2d_s *
create_UBspline_2d_s (Ugrid x_grid,   Ugrid y_grid,
		      BCtype_s   xBC, BCtype_s   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
UBspline_3d_s *
create_UBspline_3d_s (Ugrid x_grid,   Ugrid y_grid,   Ugrid z_grid,
		      BCtype_s  xBC,  BCtype_s   yBC, BCtype_s   zBC,
		      float *data);


/////////////////////////////////////
// Uniform, double precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
UBspline_1d_d *
create_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, double *data);

// Create 2D uniform single-precision, real Bspline
UBspline_2d_d *
create_UBspline_2d_d (Ugrid x_grid,   Ugrid y_grid,
		      BCtype_d   xBC, BCtype_d   yBC,
		      double *data);

// Create 3D uniform single-precision, real Bspline
UBspline_3d_d *
create_UBspline_3d_d (Ugrid x_grid,   Ugrid   y_grid,   Ugrid z_grid,
		      BCtype_d  xBC,  BCtype_d   yBC, BCtype_d   zBC,
		      double *data);


///////////////////////////////////////
// Uniform, single precision, complex//
///////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
UBspline_1d_c *
create_UBspline_1d_c (Ugrid x_grid, BCtype_c xBC, float *data);

// Create 2D uniform single-precision, real Bspline
UBspline_2d_c *
create_UBspline_2d_c (Ugrid   x_grid, Ugrid   y_grid,
		      BCtype_c   xBC, BCtype_c   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
UBspline_3d_c *
create_UBspline_3d_c (Ugrid  x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype_c  xBC, BCtype_c yBC, BCtype_c zBC,
		      float *data);

 
///////////////////////////////////////
// Uniform, double precision, complex//
///////////////////////////////////////
// Create 1D uniform double-precision, complex Bspline
UBspline_1d_z *
create_UBspline_1d_z (Ugrid x_grid, BCtype_z xBC, float *data);

// Create 2D uniform double-precision, complex Bspline
UBspline_2d_z *
create_UBspline_2d_z (Ugrid x_grid, Ugrid y_grid,
		      BCtype_z   xBC, BCtype_z   yBC,
		      float *data);

// Create 3D uniform double-precision, complex Bspline
UBspline_3d_z *
create_UBspline_3d_z (Ugrid  x_grid, Ugrid   y_grid, Ugrid z_grid,
		      BCtype_z  xBC, BCtype_z   yBC, BCtype_z zBC,
		      float *data);


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////            Spline evaluation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Value only
inline void
eval_UBspline_1d_s     (UBspline_1d_s * restrict spline, double x, 
			float* restrict val);
// Value and gradient
inline void
eval_UBspline_1d_s_vg  (UBspline_1d_s* restrict spline, double x, 
			float* restrict val, float* restrict grad);
// Value, gradient, and Laplacian
inline void
eval_UBspline_1d_s_vgl (UBspline_1d_s* restrict spline, double x, 
			float* restrict val, float* restrict grad, float* restrict lapl);
// Value, gradient, and Hessian
inline void
eval_UBspline_1d_s_vgh (UBspline_1d_s* restrict spline, double x, 
			float *val, float *grad, float *hess);

inline void
eval_UBspline_2d_s     (UBspline_2d_s* restrict spline, 
			double x, double y,
			float* restrict val);
inline void
eval_UBspline_2d_s_vg  (UBspline_2d_s* restrict spline, 
			double x, double y, 
		        float* restrict val, float* restrict grad);
inline void
eval_UBspline_2d_s_vgl (UBspline_2d_s* restrict spline, 
			double x, double y,
			float* restrict val, float* restrict grad, 
			float* restrict lapl);
inline void
eval_UBspline_2d_s_vgh (UBspline_2d_s* restrict spline, 
			double x, double y,
			float* restrict val, float* restrict grad, 
			float* restrict hess);

inline void
eval_UBspline_3d_s     (UBspline_3d_s* restrict spline, 
			double x, double y, double z,
			float* restrict val);
inline void
eval_UBspline_3d_s_vg  (UBspline_3d_s* restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad);
inline void
eval_UBspline_3d_s_vgl (UBspline_3d_s* restrict spline,
			double x, double y, double z,
			float* restrict val, float* restrict grad, 
			float* restrict lapl);
inline void
eval_UBspline_3d_s_vgh (UBspline_3d_s* restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad, 
			float* restrict hess);

// Similarly for the rest of the types.

void
destroy_Bspline (Bspline *ptr);

#endif BSPLINE_H
