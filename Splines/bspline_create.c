#include "bspline.h"

/* double A[4][4] = { 1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0, */
/* 		   1.0/2.0, 0.0/1.0, 1.0/2.0, 0.0, */
/* 		   1.0/1.0,-2.0/1.0, 1.0/1.0, 0.0 }; */


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

void
find_coefs_1d (Ugrid grid, BCtype_s bc, 
	       float *data,  int dstride,
	       float *coefs, int cstride)
{
  int M = grid.num;
  float bands[M+2][4];
  if (bc.lCode == PERIODIC) {
  }
  else {
    // Setup boundary conditions
    float abcd_left[4], abcd_right[4];
    // Left boundary
    if (bc.lCode == FLAT || bc.lCode == NATURAL)
      bc.lVal = 0.0;
    if (bc.lCode == FLAT || bc.lCode == DERIV1) {
      abcd_left[0] = -0.5     * grid.delta_inv;
      abcd_left[1] =  0.0     * grid.delta_inv; 
      abcd_left[2] =  0.5     * grid.delta_inv;
      abcd_left[3] =  bc.lVal * grid.delta_inv;
    }
    if (bc.lCode == NATURAL || bc.rCode == DERIV2) {
      abcd_left[0] = 1.0     * grid.delta_inv * grid.delta_inv;
      abcd_left[1] =-2.0     * grid.delta_inv * grid.delta_inv;
      abcd_left[2] = 1.0     * grid.delta_inv * grid.delta_inv;
      abcd_left[3] = bc.lVal * grid.delta_inv * grid.delta_inv;
    }
    
    // Right boundary
    if (bc.rCode == FLAT || bc.rCode == NATURAL)
      bc.rVal = 0.0;
    if (bc.rCode == FLAT) {
      abcd_right[0] = -0.5 * grid.delta_inv;
      abcd_right[1] =  0.0 * grid.delta_inv; 
      abcd_right[2] =  0.5 * grid.delta_inv;
      abcd_right[3] =  bc.rVal * grid.delta_inv;
    }
    if (bc.rCode == NATURAL) {
      abcd_right[0] = 1.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[1] =-2.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[2] = 1.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[3] = 
	bc.rVal * grid.delta_inv * grid.delta_inv;
    }
    for (int i=0; i<4; i++) {
      bands[0][i]   = abcd_left[i];
      bands[M+1][i] = abcd_right[i];
    }
    float basis[4] = {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0};
    for (int i=0; i<M; i++) {
      for (int j=0; j<3; j++)
	bands[i+1][j] = basis[j];
      bands[i+1][3] = data[i*dstride];
    }   
    // Now, solve for coefficients
    solve_deriv_interp_1d_s (bands, coefs, M, cstride);
  }
}

	       

Bspline*
create_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, float *data)
{
  // Create new spline
  UBspline_1d_s* restrict spline = malloc (sizeof(UBspline_1d_s));
  // Setup internal variables
  x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
  x_grid.delta_inv = 1.0/spline->x_grid.delta;
  spline->x_grid   = x_grid;

  int M = x_grid.num;
  spline->coefs = malloc (sizeof(float)*M);

  find_coefs_1d (spline->x_grid, xBC, data, 1, spline->coefs, 1);
    

  return (Bspline*) spline;
}
