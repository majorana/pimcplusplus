#ifndef BSPLINE_H
#define BSPLINE_H

// Conventions:
// Postfixes:  
// s:  single precision real
// d:  double precision real
// c:  single precision complex
// z:  double precision complex

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
////              Spline creation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/////////////////////////////////////
// Uniform, single precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_f (Ugrid x_grid, BCtype_s xBC, float *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_f (Ugrid x_grid,   Ugrid y_grid,
		      BCtype_s   xBC, BCtype_s   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_f (Ugrid x_grid,   Ugrid y_grid,   Ugrid z_grid,
		      BCtype_s  xBC,  BCtype_s   yBC, BCtype_s   zBC,
		      float *data);


/////////////////////////////////////
// Uniform, doulbe precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, double *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_d (Ugrid x_grid,   Ugrid y_grid,
		      BCtype_d   xBC, BCtype_d   yBC,
		      double *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_d (Ugrid x_grid,   Ugrid   y_grid,   Ugrid z_grid,
		      BCtype_d  xBC,  BCtype_d   yBC, BCtype_d   zBC,
		      double *data);


///////////////////////////////////////
// Uniform, single precision, complex//
///////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_c (Ugrid x_grid, BCtype_c xBC, float *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_c (Ugrid   x_grid, Ugrid   y_grid,
		      BCtype_c   xBC, BCtype_c   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_c (Ugrid  x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype_c  xBC, BCtype_c yBC, BCtype_c zBC,
		      float *data);

 
///////////////////////////////////////
// Uniform, double precision, complex//
///////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_z (Ugrid x_grid, BCtype xBC_z, float *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_z (Ugrid x_grid, Ugrid y_grid,
		      BCtype_z   xBC, BCtype_z   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
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
eval_UBspline_1d_f     (Bspline *spline, double x, float *val);
// Value and gradient
inline void
eval_UBspline_1d_f_vg  (Bspline *spline, double x, 
		       float *val, float *grad);
// Value, gradient, and Laplacian
inline void
eval_UBspline_1d_f_vgl (Bspline *spline, double x, 
			float *val, float *grad, float *lapl);
// Value, gradient, and Hessian
inline void
eval_UBspline_1d_f_vgh (Bspline *spline, double x, 
			float *val, float *grad, float *hess);

inline void
eval_UBspline_2d_f     (Bspline *spline, double x, float *val);
inline void
eval_UBspline_2d_f_vg  (Bspline *spline, double x, 
		        float *val, float *grad);
inline void
eval_UBspline_2d_f_vgl (Bspline *spline, double x, 
			float *val, float *grad, float *lapl);
inline void
eval_UBspline_2d_f_vgh (Bspline *spline, double x, 
			float *val, float *grad, float *hess);

inline void
eval_UBspline_3d_f (Bspline *spline, double x, float *val);
inline void
eval_UBspline_3d_f_vg (Bspline *spline, double x, 
		       float *val, float *grad);
inline void
eval_UBspline_3d_f_vgl (Bspline *spline, double x, 
			float *val, float *grad, float *lapl);
inline void
eval_UBspline_3d_f_vgh (Bspline *spline, double x, 
			float *val, float *grad, float *hess);

// Similarly for the rest of the types.

void
destroy_Bspline (Bspline *ptr);


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////           Bspline structure definitions            ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

///////////////////////////
// Single precision real //
///////////////////////////
typedef struct
{
  spline_code code;
  float *coefs;
  Ugrid x_grid;
  BCtype xBC;
#ifdef __SSE3__
  __m128 tx;
#else
  float tx[4];
#endif
} UBspline_1d_f;

typedef struct
{
  spline_code code;
  float *coefs;
  Ugrid x_grid, y_grid;
  BCtype xBC, yBC;
#ifdef __SSE3__
  __m128 tx, ty;
#else
  float tx[4], ty[4];
#endif
} UBspline_2d_f;

typedef struct
{
  spline_code code;
  float *coefs;
  Ugrid x_grid, y_grid, zgrid;
  BCtype xBC, yBC, zBC;
#ifdef __SSE3__
  __m128 tx, ty, tz;
#else
  float tx[4], ty[4], tz[4];
#endif
} UBspline_3d_f;


///////////////////////////
// Double precision real //
///////////////////////////
typedef struct
{
  spline_code code;
  double *coefs;
  Ugrid x_grid;
  BCtype xBC;
} UBspline_1d_d;

typedef struct
{
  spline_code code;
  double *coefs;
  Ugrid x_grid, y_grid;
  BCtype xBC, yBC;
} UBspline_2d_d;

typedef struct
{
  spline_code code;
  double *coefs;
  Ugrid x_grid, y_grid, zgrid;
  BCtype xBC, yBC, zBC;
} UBspline_3d_d;




#endif BSPLINE_H
