#include "BicubicSpline.h"


void
BicubicSpline::XUpdate(int i)
{
  Grid &x = *Xgrid;
  int N = x.NumPoints;

  Array<double,1> U(N);

  if (isnan(XStartDeriv(i))) {
    d2Fdx2(0,i)=0.0; // Use "Natural" boundary conditions--ie d^2F/dx^2 = 0 
    U(0) = 0.0;
  }
  else {
    d2Fdx2(0,i) = -0.5;
    U(0) = (3.0/(x(1)-x(0)))*((F(1,i)-F(0,i))/(x(1)-x(0))-XStartDeriv(i));
  }

  d2Fdx2(Nx-1,i) = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int j=1; j<Nx-1; j++) {
    double sig = (x(j) - x(j-1)) / (x(j+1)-x(j-1));
    double p = sig *d2Fdx2(j-1,i)+2.0;
    d2Fdx2(j,i) = (sig-1.0)/p;
    U(j) = (6.0*((F(j+1,i)-F(j,i))/(x(j+1)-x(j))-
		 (F(j,i)-F(j-1,i))/(x(j)-x(j-1))) /
	    (x(j+1)-x(j-1))-sig*U(j-1))/p;   
  }

  double Un, Qn;
  if (XEndDeriv(i) > 0.99e30)
    {
      Un = 0.0;
      Qn = 0.0;
    }
  else
    {
      Qn = 0.5;
      Un = (3.0/(x(N-1)-x(N-2)))*
	(XEndDeriv(i)-(F(N-1,i)-F(N-2,i))/(x(N-1)-x(N-2)));
    }
  
  d2Fdx2(N-1,i) =
    (Un-Qn*U(N-2))/(Qn*d2Fdx2(N-2,i)+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2Fdx2(k,i) = d2Fdx2(k,i)*d2Fdx2(k+1,i) + U(k);

  XUpToDate(i) = true;
}




void
BicubicSpline::YUpdate(int i)
{
  Grid &y = *Ygrid;
  int N = y.NumPoints;

  Array<double,1> U(N);

  if (isnan(YStartDeriv(i))) {
    d2Fdy2(i,0)=0.0; // Use "Natural" boundary conditions--ie d^2F/dy^2 = 0 
    U(0) = 0.0;
  }
  else {
    d2Fdy2(i,0) = -0.5;
    U(0) = (3.0/(y(1)-y(0)))*((F(i,1)-F(i,0))/(y(1)-y(0))-YStartDeriv(i));
  }

  d2Fdy2(i,Ny-1) = 0.0;
  
  // This part translated from Numerical Recipes.
  for (int j=1; j<Ny-1; j++) {
    double sig = (y(j) - y(j-1)) / (y(j+1)-y(j-1));
    double p = sig *d2Fdy2(i,j-1)+2.0;
    d2Fdy2(i,j) = (sig-1.0)/p;
    U(j) = (6.0*((F(i,j+1)-F(i,j))/(y(j+1)-y(j))-
		 (F(i,j)-F(i,j-1))/(y(j)-y(j-1))) /
	    (y(j+1)-y(j-1))-sig*U(j-1))/p;   
  }

  double Un, Qn;
  if (YEndDeriv(i) > 0.99e30)
    {
      Un = 0.0;
      Qn = 0.0;
    }
  else
    {
      Qn = 0.5;
      Un = (3.0/(y(N-1)-y(N-2)))*
	(YEndDeriv(i)-(F(i,N-1)-F(i,N-2))/(y(N-1)-y(N-2)));
    }
  
  d2Fdy2(i,N-1) =
    (Un-Qn*U(N-2))/(Qn*d2Fdy2(i,N-2)+1.0);
  
  for (int k=N-2; k>=0; k--)
    d2Fdy2(i,k) = d2Fdy2(i,k)*d2Fdy2(i,k+1) + U(k);

  YUpToDate(i) = true;
}


void BicubicSpline::BiUpdate()
{
  // First, update X and Y splines
  for (int i=0; i<Nx; i++)
    if (!XUpToDate(i)) XUpdate(i);
  for (int i=0; i<Ny; i++)
    if (!XUpToDate(i)) YUpdate(i);
}
  


