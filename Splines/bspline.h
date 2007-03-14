#ifndef BSPLINE_H
#define BSPLINE_H

// Conventions:
// Postfixes:  
// f:  single precision real
// d:  double precision real
// c:  single precision complex
// z:  double precision complex

typedef enum { PERIODIC, DERIV1, DERIV2, FLAT, NATURAL } bc_code;
typedef enum { U1D, U2D, U3D, NU1D, NU2D, NU3D } spline_code;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal, rVal;
} BCtype;

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

typedef struct
{
  spline_code code;
  float *data;
  Ugrid x_grid;
  BCtype xBC;
} UBspline_1d;

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
create_UBspline_1d_f (Ugrid x_grid, BCtype xBC, float *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_f (Ugrid x_grid, Ugrid y_grid,
		      BCtype   xBC, BCtype   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_f (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype  xBC,  BCtype   yBC, BCtype   zBC,
		      float *data);


/////////////////////////////////////
// Uniform, doulbe precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_d (Ugrid x_grid, BCtype xBC, double *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_d (Ugrid x_grid, Ugrid y_grid,
		      BCtype   xBC, BCtype   yBC,
		      double *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype  xBC,  BCtype   yBC, BCtype   zBC,
		      double *data);


///////////////////////////////////////
// Uniform, single precision, complex//
///////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_c (Ugrid x_grid, BCtype xBC, float *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_c (Ugrid x_grid, Ugrid y_grid,
		      BCtype   xBC, BCtype   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_c (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype  xBC,  BCtype   yBC, BCtype   zBC,
		      float *data);

 
///////////////////////////////////////
// Uniform, double precision, complex//
///////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
Bspline *
create_UBspline_1d_z (Ugrid x_grid, BCtype xBC, float *data);

// Create 2D uniform single-precision, real Bspline
Bspline *
create_UBspline_2d_z (Ugrid x_grid, Ugrid y_grid,
		      BCtype   xBC, BCtype   yBC,
		      float *data);

// Create 3D uniform single-precision, real Bspline
Bspline *
create_UBspline_3d_z (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype  xBC,  BCtype   yBC, BCtype   zBC,
		      float *data);


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////            Spline evaluation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Value only
void
eval_UBspline_1d_f     (Bspline *spline, double x, float *val);
// Value and gradient
void
eval_UBspline_1d_f_vg  (Bspline *spline, double x, 
		       float *val, float *grad);
// Value, gradient, and Laplacian
void
eval_UBspline_1d_f_vgl (Bspline *spline, double x, 
			float *val, float *grad, float *lapl);
// Value, gradient, and Hessian
void
eval_UBspline_1d_f_vgh (Bspline *spline, double x, 
			float *val, float *grad, float *hess);

void
eval_UBspline_2d_f     (Bspline *spline, double x, float *val);
void
eval_UBspline_2d_f_vg  (Bspline *spline, double x, 
		        float *val, float *grad);
void
eval_UBspline_2d_f_vgl (Bspline *spline, double x, 
			float *val, float *grad, float *lapl);
void
eval_UBspline_2d_f_vgh (Bspline *spline, double x, 
			float *val, float *grad, float *hess);

void
eval_UBspline_3d_f (Bspline *spline, double x, float *val);
void
eval_UBspline_3d_f_vg (Bspline *spline, double x, 
		       float *val, float *grad);
void
eval_UBspline_3d_f_vgl (Bspline *spline, double x, 
			float *val, float *grad, float *lapl);
void
eval_UBspline_3d_f_vgh (Bspline *spline, double x, 
			float *val, float *grad, float *hess);

// Similarly for the rest of the types.

void
destroy_Bspline (Bspline *ptr);


#endif BSPLINE_H
