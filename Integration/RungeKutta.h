#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "../Splines/Grid.h"

inline double mag (double x)
{
  return (fabs(x));
}

inline double mag (const Vec2 &x)
{
  return (fabs(x[0]) + fabs(x[1]));
}

inline double mag (const Vec3 &x)
{
  return (fabs(x[0])+fabs(x[1])+fabs(x[2]));
}

template <class IntegrandClass, class T>
class RungeKutta
{
private:
  IntegrandClass &Integrand;
public:
  inline void Integrate (const Grid &grid, int startPoint, int endPoint,
			 Array<T,1> &result, bool scale=true)
  {
    T k1, k2, k3, k4;
    T y, yplus;
    double h, x;
    int direction;
    
    if (startPoint < endPoint)
      direction = 1;
    else
      direction = -1;
    
    const static double OneSixth = 1.0/6.0;
    const static double OneThird = 1.0/3.0;
    
    int i = startPoint;
    while (i!=endPoint) {
      x = grid(i);
      h = grid(i+direction) - x;
      y = result (i);
      k1 = h * Integrand(x,       y);
      yplus = y + 0.5*k1;
      k2 = h * Integrand(x+0.5*h, yplus);
      yplus = y + 0.5*k2;
      k3 = h * Integrand(x+0.5*h, yplus);
      yplus = y + k3;
      k4 = h * Integrand(x+h, yplus);
      result(i+direction) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      if (scale && mag(result(i+direction)) > 1.0e10)
	for (int j=startPoint; j!=(endPoint+direction); j+=direction)
	  result(j) *= 1.0e-10;
      i += direction;
    }
  }


  RungeKutta(IntegrandClass &integrand) : Integrand(integrand)
  {
    // Do nothing 
  }
};


#endif
