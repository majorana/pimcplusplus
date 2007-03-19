#ifndef BSPLINE_STRUCTS_STD_H
#define BSPLINE_STRUCTS_STD_H

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
  float tx[4];
} UBspline_1d_s;

typedef struct
{
  spline_code code;
  float *coefs;
  int coefs_size, x_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
  float tx[4], ty[4];
} UBspline_2d_s;

typedef struct
{
  spline_code code;
  float *coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  float tx[4], ty[4], tz[4];
  int csize;

} UBspline_3d_s;


///////////////////////////
// Double precision real //
///////////////////////////
typedef struct
{
  spline_code code;
  double *coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_d xBC;
} UBspline_1d_d;

typedef struct
{
  spline_code code;
  double *coefs;
  int coefs_size, x_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
} UBspline_2d_d;

typedef struct
{
  spline_code code;
  double *coefs;
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
  float* restrict coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_s xBC;
  float tx[4];
} UBspline_1d_c;

typedef struct
{
  spline_code code;
  float *coefs;
  int coefs_size, x_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
  float tx[4], ty[4];
} UBspline_2d_c;

typedef struct
{
  spline_code code;
  float *coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  float tx[4], ty[4], tz[4];
  int csize;

} UBspline_3d_c;


//////////////////////////////
// Double precision complex //
//////////////////////////////
typedef struct
{
  spline_code code;
  double *coefs;
  int coefs_size;
  Ugrid x_grid;
  BCtype_d xBC;
} UBspline_1d_z;

typedef struct
{
  spline_code code;
  double *coefs;
  int coefs_size, x_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
} UBspline_2d_z;

typedef struct
{
  spline_code code;
  double *coefs;
  int coefs_size, x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
} UBspline_3d_z;


#endif
