#include "../Blitz.h"

inline double LegendrePoly (int n, double x)
{
  double Pnm2 = 1.0;
  double Pnm1 = x;
  double Pn;
  
  if (n == 0)
    return (1.0);
  else if (n == 1)
    return x;
  else 
    for (int i=2; i<=n; i++) {
      Pn = ((2*i-1)*x*Pnm1-(i-1)*Pnm2)/i;
      Pnm2 = Pnm1;
      Pnm1 = Pn;
    }
  return Pn;
}


inline void LegendrePoly (double x, Array<double,1> &Pn)
{
  int nmax = Pn.size()-1;
  if (nmax >= 0)
    Pn(0) = 1.0;
  if (nmax >= 1)
    Pn(1) = x;
  for (int n=2; n<=nmax; n++)
    Pn(n) = ((2*n-1)*x*Pn(n-1)-(n-1)*Pn(n-2))/n;
}
