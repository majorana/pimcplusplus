#ifndef RUNGA_KUTTA_H
#define RUNGA_KUTTA_H

#include "../Splines/Grid.h"

template <class T>
class FirstOrderRK
{
private:
  T &Integrand;
public:
  inline void Integrate (const Grid &grid, int StartPoint, int EndPoint,
		  Array<double,1 &result)
  {
    scalar k1, k2, k3, k4;
    scalar y, yplus;
    scalar h, x;
    int direction;
    
    if (StartPoint < EndPoint)
      direction = 1;
    else
      direction = -1;
    
    const static scalar OneSixth = 1.0/6.0;
    const static scalar OneThird = 1.0/3.0;
    
    int i = StartPoint;
    while (i!=EndPoint) {
      x = grid(i);
      h = grid(i+direction) - x;
      y = Result (i);
      k1 = h * Integrand(x,       y);
      yplus = y + 0.5*k1;
      k2 = h * Integrand(x+0.5*h, yplus);
      yplus = y + 0.5*k2;
      k3 = h * Integrand(x+0.5*h, yplus);
      yplus = y + k3;
      k4 = h * Integrand(x+h, yplus);
      Result(i+direction) = 
	y + (OneSixth*(k1+k4) + OneThird*(k2+k3));
      i += direction;
    }
  }


  RK_FirstOrder(T &integrand) : Integrand(integrand)
  {
    // Do nothing 
  }
};





#endif
