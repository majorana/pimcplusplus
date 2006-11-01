#include "TricubicBspline.h"

void
TricubicBspline::SolvePeriodicInterp (Array<double,1> data, Array<double,1> p)
{
  double ratio = 0.25;
  int N = data.size();

  Array<double,1> d(N), gamma(N), mu(N);
  d = 1.5*data;
  // First, eliminate leading coefficients
  gamma (0) = ratio;
  mu(0) = ratio;
  mu(N-1) = ratio;
  gamma(N-1) = 1.0;
  for (int row=1; row <(N-1); row++) {
    double diag = 1.0- mu(row-1)*ratio;
    double diagInv = 1.0/diag;
    gamma(row) = -ratio*gamma(row-1)*diagInv;
    mu(row) = diagInv*ratio;
    d(row)  = diagInv*(d(row)-ratio*d(row-1));
    // Last row
    d(N-1) -= mu(N-1) * d(row-1);
    gamma(N-1) -= mu(N-1)*gamma(row-1);
    mu(N-1) = -mu(N-1)*mu(row-1);
  }
  // Last row:  gamma(N-1) hold diagonal element
  mu(N-1) += ratio;
  gamma(N-1) -= mu(N-1)*(mu(N-2)+gamma(N-2));
  d(N-1) -= mu(N-1) * d(N-2);
  p(N-1) = d(N-1)/gamma(N-1);
 
  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    p(row) = d(row) - mu(row)*p(row+1) - gamma(row)*p(N-1);

}

void
TricubicBspline::SolvePeriodicInterp (Array<double,3> &data)
{
  // Do X direction
  for (int iy=0; iy<Ny; iy++)
    for (int iz=0; iz<Nz; iz++) 
      SolvePeriodicInterp(data(Range(0,Nx-1), iy, iz), P(Range(-1,Nx-2),iy, iz));
  
  // Do Y direction
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) 
      SolvePeriodicInterp(P(ix,Range(0,Ny-1), iz), P(ix, Range(-1,Ny-2), iz));
  
  // Do z direction
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) 
      SolvePeriodicInterp(P(ix,iy,Range(0,Nz-1)), P(ix, iy, Range(-1,Nz-2)));
}

void
TricubicBspline::MakePeriodic()
{
  // Now, make periodic
  for (int ix=-1; ix<(Nx+2); ix++)
    for (int iy=-1; iy<(Ny+2); iy++)
      for (int iz=-1; iz<(Nz+2); iz++)
	P(ix, iy, iz) = P((ix+Nx)%Nx, (iy+Ny)%Ny, (iz+Nz)%Nz);
}


void
TricubicBspline::Init (double xi, double xf, double yi, double yf, double zi, double zf,
		       Array<double,3> &data, bool interp, bool periodic)
{
  Nx = data.extent(0);
  Ny = data.extent(1);
  Nz = data.extent(2);
  xStart=xi;  xEnd=xf;  yStart=yi;  yEnd=yf;  zStart=zi; zEnd=zf;
  Lx    = xf-xi;  Ly    = yf-yi;  Lz    = zf-zi;
  LxInv = 1.0/Lx; LyInv = 1.0/Ly; LzInv = 1.0/Lz;
  dx = Lx/(double)Nx; dxInv = 1.0/dx;
  dy = Ly/(double)Ny; dyInv = 1.0/dy;
  dz = Lz/(double)Nz; dzInv = 1.0/dz;

  Periodic = periodic;
  Interpolating = interp;

  P.resize(Range(-1,Nx+1), Range(-1,Ny+1), Range(-1,Nz+1));

  if (Periodic) {
    if (interp)
      SolvePeriodicInterp(data);
    else
      P(Range(0,Nx-1),Range(0,Ny-1),Range(0,Nz-1)) = data;
    MakePeriodic();
  }
  else {
    cerr << "Nonperiodic Tricubic B-splines not yet supported.\n";
    abort();
  }
}

TricubicBspline::TricubicBspline()
{
  A(0,0) = -1.0/6.0; A(0,1) =  3.0/6.0; A(0,2) = -3.0/6.0; A(0,3) = 1.0/6.0;
  A(1,0) =  3.0/6.0; A(1,1) = -6.0/6.0; A(1,2) =  3.0/6.0; A(1,3) = 0.0/6.0;
  A(2,0) = -3.0/6.0; A(2,1) =  0.0/6.0; A(2,2) =  3.0/6.0; A(2,3) = 0.0/6.0;
  A(3,0) =  1.0/6.0; A(3,1) =  4.0/6.0; A(3,2) =  1.0/6.0; A(3,3) = 0.0/6.0;

  dA(0,0)= 0.0; dA(0,1)= 0.0; dA(0,2)= 0.0; dA(0,3)= 0.0;
  dA(1,0)=-0.5; dA(1,1)= 1.5; dA(1,2)=-1.5; dA(1,3)= 0.5;
  dA(2,0)= 1.0; dA(2,1)=-2.0; dA(2,2)= 1.0; dA(2,3)= 0.0;
  dA(3,0)=-0.5; dA(3,1)= 0.0; dA(3,2)= 0.5; dA(3,3)= 0.0;

  d2A(0,0)= 0.0; d2A(0,1)= 0.0; d2A(0,2)= 0.0; d2A(0,3)= 0.0;
  d2A(1,0)= 0.0; d2A(1,1)= 0.0; d2A(1,2)= 0.0; d2A(1,3)= 0.0;
  d2A(2,0)=-1.0; d2A(2,1)= 3.0; d2A(2,2)=-3.0; d2A(2,3)= 1.0;
  d2A(3,0)= 1.0; d2A(3,1)=-2.0; d2A(3,2)= 1.0; d2A(3,3)= 0.0;

  d3A(0,0)= 0.0; d3A(0,1)= 0.0; d3A(0,2)= 0.0; d3A(0,3)= 0.0;
  d3A(1,0)= 0.0; d3A(1,1)= 0.0; d3A(1,2)= 0.0; d3A(1,3)= 0.0;
  d3A(2,0)= 0.0; d3A(2,1)= 0.0; d3A(1,2)= 2.0; d3A(2,3)= 0.0;
  d3A(3,0)=-1.0; d3A(3,1)= 3.0; d3A(3,2)=-3.0; d3A(3,3)= 1.0;
}



