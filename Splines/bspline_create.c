#include "bspline.h"

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void 
solve_deriv_interp_1d_s (float bands[][4], float coefs[],
			 int M, int cstride)
{
  // Solve interpolating equations
  // First and last rows are different
  bands[0][1] /= bands[0][0];
  bands[0][2] /= bands[0][0];
  bands[0][3] /= bands[0][0];
  bands[0][0] = 1.0;
  bands[1][1] -= bands[1][0]*bands[0][1];
  bands[1][2] -= bands[1][0]*bands[0][2];
  bands[1][3] -= bands[1][0]*bands[0][3];
  bands[0][0] = 0.0;
  bands[1][2] /= bands[1][1];
  bands[1][3] /= bands[1][1];
  bands[1][1] = 1.0;
  
  // Now do rows 2 through M+1
  for (int row=2; row < (M+1); row++) {
    bands[row][1] -= bands[row][0]*bands[row-1][2];
    bands[row][3] -= bands[row][0]*bands[row-1][3];
    bands[row][2] /= bands[row][1];
    bands[row][3] /= bands[row][1];
    bands[row][0] = 0.0;
    bands[row][1] = 1.0;
  }

  // Do last row
  bands[M+1][1] -= bands[M+1][0]*bands[M-1][2];
  bands[M+1][3] -= bands[M+1][0]*bands[M-1][3];
  bands[M+1][2] -= bands[M+1][1]*bands[M][2];
  bands[M+1][3] -= bands[M+1][1]*bands[M][3];
  bands[M+1][3] /= bands[M+1][2];
  bands[M+1][2] = 1.0;

  coefs[(M+1)*cstride] = bands[M+1][3];
  // Now back substitute up
  for (int row=M; row>0; row--)
    coefs[row*cstride] = bands[row][3] - bands[row][2]*coefs[cstride*(row+1)];
  
  // Finish with first row
  coefs[0] = bands[0][3] - bands[0][1]*coefs[1*cstride] - bands[0][2]*coefs[2*cstride];
}

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void 
solve_periodic_interp_1d_s (float bands[][4], float coefs[],
			    int M, int cstride)
{
  float lastCol[M];
  // Now solve:
  // First and last rows are different
  bands[0][2] /= bands[0][1];
  bands[0][0] /= bands[0][1];
  bands[0][3] /= bands[0][1];
  bands[0][1]  = 1.0;
  bands[M-1][1] -= bands[M-1][2]*bands[0][0];
  bands[M-1][3] -= bands[M-1][2]*bands[0][3];
  bands[M-1][2]  = -bands[M-1][2]*bands[0][2];
  lastCol[0] = bands[0][0];
  
  for (int row=1; row < (M-1); row++) {
    bands[row][1] -= bands[row][0] * bands[row-1][2];
    bands[row][3] -= bands[row][0] * bands[row-1][3];
    lastCol[row]   = -bands[row][0] * lastCol[row-1];
    bands[row][0] = 0.0;
    bands[row][2] /= bands[row][1];
    bands[row][3] /= bands[row][1];
    lastCol[row]  /= bands[row][1];
    bands[row][1]  = 1.0;
    if (row < (M-2)) {
      bands[M-1][3] -= bands[M-1][2]*bands[row][3];
      bands[M-1][1] -= bands[M-1][2]*lastCol[row];
      bands[M-1][2] = -bands[M-1][2]*bands[row][2];
    }
  }

  // Now do last row
  // The [2] element and [0] element are now on top of each other 
  bands[M-1][0] += bands[M-1][2];
  bands[M-1][1] -= bands[M-1][0] * (bands[M-2][2]+lastCol[M-2]);
  bands[M-1][3] -= bands[M-1][0] *  bands[M-2][3];
  bands[M-1][3] /= bands[M-1][1];
  coefs[M*cstride] = bands[M-1][3];
  for (int row=M-2; row>=0; row--) 
    coefs[(row+1)*cstride] = 
      bands[row][3] - bands[row][2]*coefs[(row+2)*cstride] - lastCol[row]*coefs[M*cstride];
  
  coefs[0*cstride] = coefs[M*cstride];
  coefs[(M+1)*cstride] = coefs[1*cstride];
  coefs[(M+2)*cstride] = coefs[2*cstride];


}


void
find_coefs_1d (Ugrid grid, BCtype_s bc, 
	       float *data,  int dstride,
	       float *coefs, int cstride)
{
  int M = grid.num;
  float basis[4] = {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0};
  if (bc.lCode == PERIODIC) {
    float bands[M][4];
    for (int i=0; i<M; i++) {
      bands[i][0] = basis[0];
      bands[i][1] = basis[1];
      bands[i][2] = basis[2];
      bands[i][3] = data[i*dstride];
    }
    solve_periodic_interp_1d_s (bands, coefs, M, cstride);
  }
  else {
    float bands[M+2][4];
    // Setup boundary conditions
    float abcd_left[4], abcd_right[4];
    // Left boundary
    if (bc.lCode == FLAT || bc.lCode == NATURAL)
      bc.lVal = 0.0;
    if (bc.lCode == FLAT || bc.lCode == DERIV1) {
      abcd_left[0] = -0.5 * grid.delta_inv;
      abcd_left[1] =  0.0 * grid.delta_inv; 
      abcd_left[2] =  0.5 * grid.delta_inv;
      abcd_left[3] =  bc.lVal;
    }
    if (bc.lCode == NATURAL || bc.lCode == DERIV2) {
      abcd_left[0] = 1.0 * grid.delta_inv * grid.delta_inv;
      abcd_left[1] =-2.0 * grid.delta_inv * grid.delta_inv;
      abcd_left[2] = 1.0 * grid.delta_inv * grid.delta_inv;
      abcd_left[3] = bc.lVal;
    }
    
    // Right boundary
    if (bc.rCode == FLAT || bc.rCode == NATURAL)
      bc.rVal = 0.0;
    if (bc.rCode == FLAT || bc.rCode == DERIV1) {
      abcd_right[0] = -0.5 * grid.delta_inv;
      abcd_right[1] =  0.0 * grid.delta_inv; 
      abcd_right[2] =  0.5 * grid.delta_inv;
      abcd_right[3] =  bc.rVal;
    }
    if (bc.rCode == NATURAL || bc.rCode == DERIV2) {
      abcd_right[0] = 1.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[1] =-2.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[2] = 1.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[3] = bc.rVal;
    }
    for (int i=0; i<4; i++) {
      bands[0][i]   = abcd_left[i];
      bands[M+1][i] = abcd_right[i];
    }
    for (int i=0; i<M; i++) {
      for (int j=0; j<3; j++)
	bands[i+1][j] = basis[j];
      bands[i+1][3] = data[i*dstride];
    }   
    // Now, solve for coefficients
    solve_deriv_interp_1d_s (bands, coefs, M, cstride);
  }
}

	       

UBspline_1d_s*
create_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, float *data)
{
  // Create new spline
  UBspline_1d_s* restrict spline = malloc (sizeof(UBspline_1d_s));
  // Setup internal variables
  int M = x_grid.num;
  int N;

  if (xBC.lCode == PERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    N = M+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    N = M+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
  spline->coefs = malloc (sizeof(float)*N);
  find_coefs_1d (spline->x_grid, xBC, data, 1, spline->coefs, 1);
    
  return spline;
}


UBspline_2d_s*
create_UBspline_2d_s (Ugrid x_grid, Ugrid y_grid,
		      BCtype_s xBC, BCtype_s yBC, float *data)
{
  // Create new spline
  UBspline_2d_s* restrict spline = malloc (sizeof(UBspline_2d_s));
  // Setup internal variables
  int Mx = x_grid.num;
  int My = y_grid.num;
  int Nx, Ny;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;
  spline->x_stride = Ny;

  spline->coefs = malloc (sizeof(float)*Nx*Ny);

  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) {
    int doffset = iy;
    int coffset = iy;
    find_coefs_1d (spline->x_grid, xBC, data+doffset, My,
		   spline->coefs+coffset, Ny);
  }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) {
    int doffset = ix*Ny;
    int coffset = ix*Ny;
    find_coefs_1d (spline->y_grid, yBC, spline->coefs+doffset, 1, 
		   spline->coefs+coffset, 1);
  }
  
  return spline;
}
