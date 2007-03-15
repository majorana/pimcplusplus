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
/* 1D double-precision, real evaulation functions           */
/************************************************************/



#endif
