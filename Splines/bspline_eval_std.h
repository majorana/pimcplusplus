#ifndef BSPLINE_EVAL_STD_H
#define BSPLINE_EVAL_STD_H

extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_s (UBspline_1d_s * restrict spline, 
		    double x, float* restrict val)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and first derivative */
inline void
eval_UBspline_1d_s_vg (UBspline_1d_s * restrict spline, double x, 
		     float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_s_vgl (UBspline_1d_s * restrict spline, double x, 
			float* restrict val, float* restrict grad,
			float* restrict lapl)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv * 
    (coefs[i+0]*(d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3])+
     coefs[i+1]*(d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3])+
     coefs[i+2]*(d2Af[10]*tp[2] + d2Af[11]*tp[3])+
     coefs[i+3]*(d2Af[14]*tp[2] + d2Af[15]*tp[3]));
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
		    double x, double y, float* restrict val)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  
  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	  a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	  a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	  a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
#undef C

}


/* Value and gradient */
inline void
eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
		       double x, double y, 
		       float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0]  = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  
  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	     a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	     a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	     a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[0] = (da[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	     da[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	     da[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	     da[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[1] = (a[0]*(C(0,0)*db[0]+C(0,1)*db[1]+C(0,2)*db[2]+C(0,3)*db[3])+
	     a[1]*(C(1,0)*db[0]+C(1,1)*db[1]+C(1,2)*db[2]+C(1,3)*db[3])+
	     a[2]*(C(2,0)*db[0]+C(2,1)*db[1]+C(2,2)*db[2]+C(2,3)*db[3])+
	     a[3]*(C(3,0)*db[0]+C(3,1)*db[1]+C(3,2)*db[2]+C(3,3)*db[3]));
#undef C

}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_s_vgl (UBspline_2d_s * restrict spline, 
			double x, double y, float* restrict val, 
			float* restrict grad, float* restrict lapl)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  
  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	     a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	     a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	     a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[0] = (da[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	     da[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	     da[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	     da[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[1] = (a[0]*(C(0,0)*db[0]+C(0,1)*db[1]+C(0,2)*db[2]+C(0,3)*db[3])+
	     a[1]*(C(1,0)*db[0]+C(1,1)*db[1]+C(1,2)*db[2]+C(1,3)*db[3])+
	     a[2]*(C(2,0)*db[0]+C(2,1)*db[1]+C(2,2)*db[2]+C(2,3)*db[3])+
	     a[3]*(C(3,0)*db[0]+C(3,1)*db[1]+C(3,2)*db[2]+C(3,3)*db[3]));
  *lapl   = ((a[0]*(C(0,0)*d2b[0]+C(0,1)*d2b[1]+C(0,2)*d2b[2]+C(0,3)*d2b[3])+
	      a[1]*(C(1,0)*d2b[0]+C(1,1)*d2b[1]+C(1,2)*d2b[2]+C(1,3)*d2b[3])+
	      a[2]*(C(2,0)*d2b[0]+C(2,1)*d2b[1]+C(2,2)*d2b[2]+C(2,3)*d2b[3])+
	      a[3]*(C(3,0)*d2b[0]+C(3,1)*d2b[1]+C(3,2)*d2b[2]+C(3,3)*d2b[3])) +
	     (d2a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	      d2a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	      d2a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	      d2a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3])));
  
#undef C

}


/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_s_vgh (UBspline_2d_s * restrict spline, 
			double x, double y, float* restrict val, 
			float* restrict grad, float* restrict hess)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  
  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    (  a[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
	       a[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
	       a[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
	       a[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  grad[0] = ( da[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
	      da[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
	      da[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
	      da[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  grad[1] = (  a[0]*(C(0,0)* db[0]+C(0,1)* db[1]+C(0,2)* db[2]+C(0,3)* db[3])+
	       a[1]*(C(1,0)* db[0]+C(1,1)* db[1]+C(1,2)* db[2]+C(1,3)* db[3])+
	       a[2]*(C(2,0)* db[0]+C(2,1)* db[1]+C(2,2)* db[2]+C(2,3)* db[3])+
	       a[3]*(C(3,0)* db[0]+C(3,1)* db[1]+C(3,2)* db[2]+C(3,3)* db[3]));
  hess[0] = (d2a[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
	     d2a[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
	     d2a[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
	     d2a[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  hess[1] = ( da[0]*(C(0,0)* db[0]+C(0,1)* db[1]+C(0,2)* db[2]+C(0,3)* db[3])+
	      da[1]*(C(1,0)* db[0]+C(1,1)* db[1]+C(1,2)* db[2]+C(1,3)* db[3])+
	      da[2]*(C(2,0)* db[0]+C(2,1)* db[1]+C(2,2)* db[2]+C(2,3)* db[3])+
	      da[3]*(C(3,0)* db[0]+C(3,1)* db[1]+C(3,2)* db[2]+C(3,3)* db[3]));
  hess[3] = (  a[0]*(C(0,0)*d2b[0]+C(0,1)*d2b[1]+C(0,2)*d2b[2]+C(0,3)*d2b[3])+
	       a[1]*(C(1,0)*d2b[0]+C(1,1)*d2b[1]+C(1,2)*d2b[2]+C(1,3)*d2b[3])+
	       a[2]*(C(2,0)*d2b[0]+C(2,1)*d2b[1]+C(2,2)*d2b[2]+C(2,3)*d2b[3])+
	       a[3]*(C(3,0)*d2b[0]+C(3,1)*d2b[1]+C(3,2)*d2b[2]+C(3,3)*d2b[3]));
  hess[2] = hess[1];
  
#undef C

}

#endif
