#include "FastCubicSpline.h"

void
FastCubicSpline::Init (double xStart, double xEnd, Array<double,1> NewYs,
		       double startDeriv, double endDeriv)
{
  int nx = NewYs.size();
  F.resize(nx);
  Xstart = xStart;
  Xend   = xEnd;
  dx = (xEnd-xStart)/(double)(nx-1);
  dxInv = 1.0/dx;
  UpdateNatural();
}

void
FastCubicSpline::UpdatePeriodic()
{

}

void
FastCubicSpline::UpdateFixed()
{

}

void
FastCubicSpline::UpdateNatural()
{
  if (Periodic) {
    UpdatePeriodic();
    return;
  }
  Array<double,1> mu(F.size());
  int Nx = F.size();
  // Set up tridiagonal set of equations
  // Initialize RHS of equations
  F(0)[1] = 1.5*(F(1)[0]-F(0)[0])*dxInv;
  F(Nx-1)[1] = 1.5*(F(Nx-1)[0]-F(Nx-2)[0])
    *dxInv;
  mu(0) = 0.5;
  
  // Solve tri-diagonal set of equations.  First eliminate lower
  // elements.
  for (int j=1; j<(Nx-1); j++) {
    double lambda = 0.25;
    mu(j) = 0.5 - lambda;
    F(j)[1] = 3.0*(lambda*(F(j)[0]-F(j-1)[0])*dxInv+ 
		   mu(j) *(F(j+1)[0]-F(j)[0])*dxInv);
    double cj = 1.0 - lambda * mu(j-1);
    F(j)[1] -= lambda * F(j-1)[1];
    mu(j) /= cj;
    F(j)[1] /= cj;
  }
  
  // Last element
  int j = Nx-1;
  double lambda = 0.5;
  mu(j) = 0.5 - lambda;
  F(j)[1] = 
    3.0*(lambda*(F(j)[0]-F(j-1)[0])*dxInv);
  double cj = 1.0 - lambda * mu(j-1);
  F(j)[1] -= lambda * F(j-1)[1];
  mu(j) /= cj;
  F(j)[1] /= cj;
  
  // Now last d/dx is correct.  We proceed upward, back substituting.
  for (j=Nx-2; j>=0; j--) {
    F(j)[1] -= mu(j) * F(j+1)[1];
  }
}

