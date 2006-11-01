#ifndef CUBIC_B_SPLINES
#define CUBIC_B_SPLINES

#include <blitz/array.h>
#include <blitz/tinymat.h>

using namespace blitz;

class CubicBspline
{
private:
  TinyMatrix<double,4,4> A;
  double GridStart, GridEnd, GridDelta, GridDeltaInv;
  bool Interpolating, Periodic;

  // The control points
  Array<double,1> P;

  // Interpolating solvers:
  void SolvePeriodicInterp(Array<double,1> &data);  

public:
  inline double GetControlPoint (int i);
  void Set(double start, double end, Array<double,1> &data, 
	   bool interpolating=true, bool periodic=true);
  inline double operator()(double x) const;
  inline double Deriv     (double x) const;
  inline double Deriv2    (double x) const;
  inline double Deriv3    (double x) const;

  CubicBspline();
};

inline double 
CubicBspline::operator()(double x) const
{
  double fi = (x-GridStart)*GridDeltaInv;
  double ipart, t;
  t = modf (fi, &ipart);
  int i, im1, ip1, ip2;
  i = (int) ipart;
  double tv[4];
  if (Periodic) {
    int N = P.size();
    i   = i % P.size();
    im1 = (i-1+N)%N;
    ip1 = (i+1)%N;
    ip2 = (i+2)%N;
  }
  else {
    i += 1;
    im1 = i-1;
    ip1 = i+1;
    ip2 = i+2;
  }
  tv[0] = t*t*t;
  tv[1] = t*t;
  tv[2] = t;
  tv[3] = 1.0;
//   return (P(im1)*(A(0,0)*tv[0]+A(1,0)*tv[1]+A(2,0)*tv[2]+A(3,0)*tv[3])+
// 	  P(i  )*(A(0,1)*tv[0]+A(1,1)*tv[1]+A(2,1)*tv[2]+A(3,1)*tv[3])+
// 	  P(ip1)*(A(0,2)*tv[0]+A(1,2)*tv[1]+A(2,2)*tv[2]+A(3,2)*tv[3])+
// 	  P(ip2)*(A(0,3)*tv[0]+A(1,3)*tv[1]+A(2,3)*tv[2]+A(3,3)*tv[3]));
  return (tv[0]*(A(0,0)*P(im1)+A(0,1)*P(i)+A(0,2)*P(ip1)+A(0,3)*P(ip2))+
	  tv[1]*(A(1,0)*P(im1)+A(1,1)*P(i)+A(1,2)*P(ip1)+A(1,3)*P(ip2))+
	  tv[2]*(A(2,0)*P(im1)+A(2,1)*P(i)+A(2,2)*P(ip1)+A(2,3)*P(ip2))+
	  tv[3]*(A(3,0)*P(im1)+A(3,1)*P(i)+A(3,2)*P(ip1)+A(3,3)*P(ip2)));
}

#endif
