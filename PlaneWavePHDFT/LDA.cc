#include "ConjGradMPI.h"
#include "../MatrixOps/MatrixOps.h"
#include "../MPI/Communication.h"
#include "../DFT/Functionals.h"

void
ConjGradMPI::SolveLDA()
{
  int numSCIters = 0;
  bool SC = false;
  while (!SC) {
    CalcOccupancies();
    CalcChargeDensity();
    MixChargeDensity();
    CalcVHXC();
    SetTolerance(1.0e-6);

    numSCIters ++;
  }

}

void
ConjGradMPI::CalcOccupancies()
{


}

void
ConjGradMPI::CalcChargeDensity()
{
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  if (TempRho.shape() != TinyVector<int,3>(nx,ny,nz))
    TempRho.resize(nx,ny,nz);
  if (NewRho.shape() != TinyVector<int,3>(nx,ny,nz))
    NewRho.resize(nx,ny,nz);
  TempRho = 0.0;
  zVec band;

  /// We assume that each of the bands is already normalized.
  for (int bi=MyFirstBand; bi<MyLastBand; bi++) {
    band.reference(Bands(bi, Range::all()));
    FFT.PutkVec (band);
    FFT.k2r();
    for (int ix=0; ix<nx; ix++)
      for (int iy=0; iy<ny; iy++)
	for (int iz=0; iz<nz; iz++)
	  TempRho(ix,iy,iz) += Occupancies(bi)*norm (FFT.rBox(ix,iy,iz));
  }
  /// First, sum over all the band procs
  BandComm.AllSum (TempRho, NewRho);
  TempRho = NewRho;
  // Now sum over all the kPoints
  if (BandComm.MyProc() == 0) 
    kComm.AllSum(TempRho, NewRho);
  // And broadcast;
  BandComm.Broadcast(0, NewRho);
  /// Multiply by the right prefactor;
  double volInv = 1.0/H.GVecs.GetBoxVol();
  double prefactor = 2.0*volInv/((double)kComm.NumProcs()*(nx*ny*nz));
  NewRho *= prefactor;
}



/// Calculate the hartree, exchange, and correlation potentials.
/// We assume that upon entry, Rho_r contains the charge density 
/// in real space and Rho_G contains it in reciporical space.
/// On exit, VH will contain the Hartree potential and VXC will
/// contain the exchange-correlation potential, both in real space.
void 
ConjGradMPI::CalcVHXC()
{
  ///////////////////////
  // Hartree potential //
  ///////////////////////
  h_G = Rho_G;
  // Compute V_H in reciporical space
  double volInv = 1.0/H.GVecs.GetBoxVol();
  double hartreeTerm = 0.0;
  double prefact = (4.0*M_PI/H.GVecs.GetBoxVol());
  EH = 0.0;
  for (int i=0; i<h_G.size(); i++) {
    h_G(i) *= prefact*H.GVecs.DeltaGInv2(i);
    EH += norm(h_G(i)) * H.GVecs.DeltaGInv2(i);
  }
  EH *= 0.5*prefact;

  // FFT back to real space
  FFT.PutkVec(h_G);
  FFT.k2r();
  copy (real(FFT.rBox), VH);

  ////////////////////////////////////
  // Exchange-correlation potential //
  ////////////////////////////////////
  EXC = 0.0;
  double exc;
  for (int ix=0; ix<Rho_r.extent(0); ix++)
    for (int iy=0; iy<Rho_r.extent(1); iy++)
      for (int iz=0; iz<Rho_r.extent(2); iz++) {
	FortranExCorr (Rho_r(ix,iy,iz), exc, VXC(ix,iy,iz));
	EXC += exc;
      }
}

// double
// ConjGradMPI::CalcHartreeTerm(int band)
// {
//   zVec psi;
//   psi.reference (Bands(band, Range::all()));
//   FFT.PutkVec (Phip);
//   FFT.k2r();
//   copy (FFT.rBox, Phip_r);
// //   for (int ix=0; ix<Nx; ix++)
// //     for (int iy=0; iy<Ny; iy++)
// //       for (int iz=0; iz<Nz; iz++)
// // 	Phip_r(ix,iy,iz) = FFT.rBox(ix,iy,iz);
//   FFT.PutkVec (psi);
//   FFT.k2r();
  
//   for (int ix=0; ix<Nx; ix++)
//     for (int iy=0; iy<Ny; iy++)
//       for (int iz=0; iz<Nz; iz++)
// 	FFT.rBox(ix,iy,iz) *= conj(Phip_r(ix,iy,iz));
//   FFT.r2k();
//   FFT.GetkVec(h_G);
//   double volInv = 1.0/H.GVecs.GetBoxVol();
//   double e2_over_eps0 = 4.0*M_PI; // ???????
//   double hartreeTerm = 0.0;
//   for (int i=0; i<h_G.size(); i++)
//     hartreeTerm += norm(h_G(i))*H.GVecs.DeltaGInv2(i);
//   hartreeTerm *= 2.0 * e2_over_eps0 * volInv;
//   return hartreeTerm;
// }


