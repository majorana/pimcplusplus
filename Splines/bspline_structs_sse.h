#ifndef BSPLINE_STRUCTS_SSE_H
#define BSPLINE_STRUCTS_SSE_H

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>

#ifdef __cplusplus
typedef complex<float>  complex_float;
typedef complex<double> complex_double;
#else
#include <complex.h>
typedef complex float  complex_float;
typedef complex double complex_double;
#endif

///////////////////////////
// Single precision real //
///////////////////////////
typedef struct
{
  spline_code code;
  float* restrict coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_s xBC;
  __m128 tx;
} UBspline_1d_s;

typedef struct
{
  spline_code code;
  float* restrict coefs;
  int x_stride;
  Ugrid coefs_size, x_grid, y_grid;
  BCtype_s xBC, yBC;
  __m128 tx, ty;
} UBspline_2d_s;

typedef struct
{
  spline_code code;
  float* restrict coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  __m128 _tpx, _tpy, _tpz;
  __m128 _delta_inv;
} UBspline_3d_s;


///////////////////////////
// Double precision real //
///////////////////////////
typedef struct
{
  spline_code code;
  double* restrict coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_d xBC;
} UBspline_1d_d;

typedef struct
{
  spline_code code;
  double* restrict coefs;
  int coefs_size, x_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
} UBspline_2d_d;

typedef struct
{
  spline_code code;
  double* restrict coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
} UBspline_3d_d;


//////////////////////////////
// Single precision complex //
//////////////////////////////
typedef struct
{
  spline_code code;
  complex_float* restrict coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_s xBC;
  __m128 tx;
} UBspline_1d_c;

typedef struct
{
  spline_code code;
  complex_float* restrict coefs;
  int x_stride;
  Ugrid coefs_size, x_grid, y_grid;
  BCtype_s xBC, yBC;
  __m128 tx, ty;
} UBspline_2d_c;

typedef struct
{
  spline_code code;
  complex_float* restrict coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  __m128 _tpx, _tpy, _tpz;
  __m128 _delta_inv;
} UBspline_3d_c;


//////////////////////////////
// Double precision complex //
//////////////////////////////
typedef struct
{
  spline_code code;
  complex_double* restrict coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_d xBC;
} UBspline_1d_z;

typedef struct
{
  spline_code code;
  complex_double* restrict coefs;
  int coefs_size, x_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
} UBspline_2d_z;

typedef struct
{
  spline_code code;
  complex_double* restrict coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
} UBspline_3d_z;



#endif
