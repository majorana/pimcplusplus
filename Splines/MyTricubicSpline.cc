#include "MyTricubicSpline.h"

void MyTricubicSpline::Update()
{
  // Do dF/dx
  UpdateX(0, 1);
  // Do dF/dy
  UpdateY(0, 2);
  // Do dF/dy
  UpdateZ(0, 3);
  // Do d2F/dxdy
  UpdateY(1, 4);
  // Do d2F/dxdz
  UpdateZ(1, 5);
  // Do d2F/dydz
  UpdateZ(2, 6);
  // Do d3F/dxdydz
  UpdateZ(4, 7);
  UpToDate=true;
}

void MyTricubicSpline::UpdateX(int source, int dest)
{
  Grid &x = *Xgrid;
  Array<double,1> mu(Nx);
  
  // Loop over all y and z
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(0,iy,iz)[dest] = 1.5*(F(1,iy,iz)[source]-F(0,iy,iz)[source])/(x(1)-x(0));
      F(Nx-1,iy,iz)[dest] = 1.5*(F(Nx-1,iy,iz)[source]-F(Nx-2,iy,iz)[source])
	/(x(Nx-1)-x(Nx-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nx-1); j++) {
	double lambda = 0.5*(x(j+1)-x(j))/(x(j+1)-x(j-1));
	mu(j) = 0.5 - lambda;
	F(j,iy,iz)[dest] = 
	  3.0*(lambda*(F(j,iy,iz)[source]-F(j-1,iy,iz)[source])/(x(j)-x(j-1))+ 
	       mu(j) *(F(j+1,iy,iz)[source]-F(j,iy,iz)[source])/(x(j+1)-x(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(j,iy,iz)[dest] -= lambda * F(j-1,iy,iz)[dest];
	mu(j) /= cj;
	F(j,iy,iz)[dest] /= cj;
      }
      
      // Last element
      int j = Nx-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(j,iy,iz)[dest] = 
	3.0*(lambda*(F(j,iy,iz)[source]-F(j-1,iy,iz)[source])/(x(j)-x(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz)[dest] -= lambda * F(j-1,iy,iz)[dest];
      mu(j) /= cj;
      F(j,iy,iz)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nx-2; j>=0; j--) {
	F(j,iy,iz)[dest] -= mu(j) * F(j+1,iy,iz)[dest];
      }
    }      
}


void MyTricubicSpline::UpdateY(int source, int dest)
{
  Grid &y = *Ygrid;
  Array<double,1> mu(Ny);
  
  // Loop over all y and z
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,0,iz)[dest] = 1.5*(F(ix,1,iz)[source]-F(ix,0,iz)[source])/(y(1)-y(0));
      F(ix,Ny-1,iz)[dest] = 1.5*(F(ix,Ny-1,iz)[source]-F(ix,Ny-2,iz)[source])
	/(y(Ny-1)-y(Ny-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Ny-1); j++) {
	double lambda = 0.5*(y(j+1)-y(j))/(y(j+1)-y(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,j,iz)[dest] = 
	  3.0*(lambda*(F(ix,j,iz)[source]-F(ix,j-1,iz)[source])/(y(j)-y(j-1))+ 
	       mu(j) *(F(ix,j+1,iz)[source]-F(ix,j,iz)[source])/(y(j+1)-y(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,j,iz)[dest] -= lambda * F(ix,j-1,iz)[dest];
	mu(j) /= cj;
	F(ix,j,iz)[dest] /= cj;
      }
      
      // Last element
      int j = Ny-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,j,iz)[dest] = 
	3.0*(lambda*(F(ix,j,iz)[source]-F(ix,j-1,iz)[source])/(y(j)-y(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,j,iz)[dest] -= lambda * F(ix,j-1,iz)[dest];
      mu(j) /= cj;
      F(ix,j,iz)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Ny-2; j>=0; j--) {
	F(ix,j,iz)[dest] -= mu(j) * F(ix,j+1,iz)[dest];
      }
    }      
}


void MyTricubicSpline::UpdateZ(int source, int dest)
{
  Grid &z = *Zgrid;
  Array<double,1> mu(Nz);
  
  // Loop over all y and z
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) {
      // Set up tridiagonal set of equations
      // Initialize RHS of equations
      F(ix,iy,0)[dest] = 1.5*(F(ix,iy,1)[source]-F(ix,iy,0)[source])/(z(1)-z(0));
      F(ix,iy,Nz-1)[dest] = 1.5*(F(ix,iy,Nz-1)[source]-F(ix,iy,Nz-2)[source])
	/(z(Nz-1)-z(Nz-2));
      mu(0) = 0.5;

      // Solve tri-diagonal set of equations.  First eliminate lower
      // elements.
      for (int j=1; j<(Nz-1); j++) {
	double lambda = 0.5*(z(j+1)-z(j))/(z(j+1)-z(j-1));
	mu(j) = 0.5 - lambda;
	F(ix,iy,j)[dest] = 
	  3.0*(lambda*(F(ix,iy,j)[source]-F(ix,iy,j-1)[source])/(z(j)-z(j-1))+ 
	       mu(j) *(F(ix,iy,j+1)[source]-F(ix,iy,j)[source])/(z(j+1)-z(j)));
	double cj = 1.0 - lambda * mu(j-1);
	F(ix,iy,j)[dest] -= lambda * F(ix,iy,j-1)[dest];
	mu(j) /= cj;
	F(ix,iy,j)[dest] /= cj;
      }
      
      // Last element
      int j = Nz-1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,iy,j)[dest] = 
	3.0*(lambda*(F(ix,iy,j)[source]-F(ix,iy,j-1)[source])/(z(j)-z(j-1)));
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,iy,j)[dest] -= lambda * F(ix,iy,j-1)[dest];
      mu(j) /= cj;
      F(ix,iy,j)[dest] /= cj;
      
      // Now last d/dx is correct.  We proceed upward, back substituting.
      for (j=Nz-2; j>=0; j--) {
	F(ix,iy,j)[dest] -= mu(j) * F(ix,iy,j+1)[dest];
      }
    }      
}


