#include "MultiTricubicSpline3.h"

void MultiTricubicSpline::Update()
{
  for (int i=0; i<N; i++) {
    // Do dF/dx
    UpdateX(0, 1, i);
    // Do dF/dy
    UpdateY(0, 2, i);
    // Do dF/dy
    UpdateZ(0, 3, i);
    // Do d2F/dxdy
    UpdateY(1, 4, i);
    // Do d2F/dxdz
    UpdateZ(1, 5, i);
    // Do d2F/dydz
    UpdateZ(2, 6, i);
    // Do d3F/dxdydz
    UpdateZ(4, 7, i);
  }
  UpToDate=true;
}

void MultiTricubicSpline::UpdateX(int source, int dest, int i)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  
  // Loop over all y and z
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(0,iy,iz,i)[dest] = 1.5*(F(1,iy,iz,i)[source]-F(0,iy,iz,i)[source])/(x(1)-x(0));
      F(Nx-1,iy,iz,i)[dest] = 1.5*(F(Nx-1,iy,iz,i)[source]-F(Nx-2,iy,iz,i)[source])
	/(x(Nx-1)-x(Nx-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nx-1); j++) {
	double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
	mu(j) = 0.5 - lambda;
	F(j,iy,iz,i)[dest] = 
	  3.0*(lambda*(F(j,iy,iz,i)[source]-F(j-1,iy,iz,i)[source])/(x(j)-x(j-1))+ 
	       mu(j) *(F(j+1,iy,iz,i)[source]-F(j,iy,iz,i)[source])/(x(j+1)-x(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(j,iy,iz,i)[dest] -= lambda * F(j-1,iy,iz,i)[dest];
	mu(j) /= cj;
	F(j,iy,iz,i)[dest] /= cj;
      }
      
      // Last element
      int j = Nx-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(j,iy,iz,i)[dest] = 
	3.0*(lambda*(F(j,iy,iz,i)[source]-F(j-1,iy,iz,i)[source])/(x(j)-x(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz,i)[dest] -= lambda * F(j-1,iy,iz,i)[dest];
      mu(j) /= cj;
      F(j,iy,iz,i)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nx-2; j>=0; j--) {
	F(j,iy,iz,i)[dest] -= mu(j) * F(j+1,iy,iz,i)[dest];
      }
    }      
}


void MultiTricubicSpline::UpdateY(int source, int dest, int i)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);
  
  // Loop over all x and z
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,0,iz,i)[dest] = 1.5*(F(ix,1,iz,i)[source]-F(ix,0,iz,i)[source])/(y(1)-y(0));
      F(ix,Ny-1,iz,i)[dest] = 1.5*(F(ix,Ny-1,iz,i)[source]-F(ix,Ny-2,iz,i)[source])
	/(y(Ny-1)-y(Ny-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Ny-1); j++) {
	double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,j,iz,i)[dest] = 
	  3.0*(lambda*(F(ix,j,iz,i)[source]-F(ix,j-1,iz,i)[source])/(y(j)-y(j-1))+ 
	       mu(j) *(F(ix,j+1,iz,i)[source]-F(ix,j,iz,i)[source])/(y(j+1)-y(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,j,iz,i)[dest] -= lambda * F(ix,j-1,iz,i)[dest];
	mu(j) /= cj;
	F(ix,j,iz,i)[dest] /= cj;
      }
      
      // Last element
      int j = Ny-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,j,iz,i)[dest] = 
	3.0*(lambda*(F(ix,j,iz,i)[source]-F(ix,j-1,iz,i)[source])/(y(j)-y(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,j,iz,i)[dest] -= lambda * F(ix,j-1,iz,i)[dest];
      mu(j) /= cj;
      F(ix,j,iz,i)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Ny-2; j>=0; j--) {
	F(ix,j,iz,i)[dest] -= mu(j) * F(ix,j+1,iz,i)[dest];
      }
    }      
}


void MultiTricubicSpline::UpdateZ(int source, int dest, int i)
{
  Grid &z = *Zgrid;
  Array<double,1> mu(Nz);
  
  // Loop over all x and y
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,iy,0,i)[dest] = 1.5*(F(ix,iy,1,i)[source]-F(ix,iy,0,i)[source])/(z(1)-z(0));
      F(ix,iy,Nz-1,i)[dest] = 1.5*(F(ix,iy,Nz-1,i)[source]-F(ix,iy,Nz-2,i)[source])
	/(z(Nz-1)-z(Nz-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nz-1); j++) {
	double lambda = 0.5*(z(j+1)-z(j))/(z(j+1)-z(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,iy,j,i)[dest] = 
	  3.0*(lambda*(F(ix,iy,j,i)[source]-F(ix,iy,j-1,i)[source])/(z(j)-z(j-1))+ 
	       mu(j) *(F(ix,iy,j+1,i)[source]-F(ix,iy,j,i)[source])/(z(j+1)-z(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,iy,j,i)[dest] -= lambda * F(ix,iy,j-1,i)[dest];
	mu(j) /= cj;
	F(ix,iy,j,i)[dest] /= cj;
      }
      
      // Last element
      int j = Nz-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,iy,j,i)[dest] = 
	3.0*(lambda*(F(ix,iy,j,i)[source]-F(ix,iy,j-1,i)[source])/(z(j)-z(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,iy,j,i)[dest] -= lambda * F(ix,iy,j-1,i)[dest];
      mu(j) /= cj;
      F(ix,iy,j,i)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nz-2; j>=0; j--) {
	F(ix,iy,j,i)[dest] -= mu(j) * F(ix,iy,j+1,i)[dest];
      }
    }      
}


