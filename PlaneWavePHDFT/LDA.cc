/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "PlaneWavesMPI.h"
#include "../MatrixOps/MatrixOps.h"
#include "../MPI/Communication.h"
#include "../DFT/Functionals.h"

void
MPISystemClass::InitLDA()
{
  ChargeMixer = new KerkerMixerClass (FFT);
  int NDelta = GVecs.DeltaSize();
  FFT.GetDims    (Nx, Ny, Nz);
  NewRho.resize  (Nx, Ny, Nz);
  Rho_r1.resize  (Nx, Ny, Nz);
  Rho_r2.resize  (Nx, Ny, Nz);
  TempRho.resize (Nx, Ny, Nz);
  VH.resize      (Nx, Ny, Nz);
  VXC.resize     (Nx, Ny, Nz);
  VHXC.resize    (Nx, Ny, Nz);
  Rho_r.resize   (Nx, Ny, Nz);
  h_G.resize     (NDelta);
  Rho_G.resize   (NDelta);
  Hmat.resize(NumBands, NumBands);
  EigVecs.resize(NumBands, NumBands);
  EigVals.resize(NumBands);
  RotBands.resize(NumBands, GVecs.size());
  Occupancies.resize(NumBands);
  CalcRadialChargeDensity();
  Smearer.SetOrder(2);
  Smearer.SetWidth(0.02);
  if (BandComm.MyProc() == 0) {
    Numk = kComm.NumProcs();
    Myk  = kComm.MyProc();
  }
  BandComm.Broadcast (0, Numk);
  BandComm.Broadcast (0, Myk);
}

#include "Hamiltonians.h"


void
MPISystemClass::SubspaceRotate()
{
  //  ApplyH();
  zVec c, Hc;
  for (int bi=0; bi<Bands.rows(); bi++) {
    c.reference (Bands(bi,Range::all()));
    for (int bj=0; bj<Bands.rows(); bj++) {
      Hc.reference (HBands(bj, Range::all()));
      Hmat(bi,bj) = conjdot (c, Hc);
    }
  }
  SymmEigenPairs(Hmat, Hmat.rows(), EigVals, EigVecs);
  // Now, compute the new rotated bands.
  RotBands = 0.0;
  for (int bi=0; bi<NumBands; bi++)
    for (int bj=0; bj<NumBands; bj++)
      RotBands(bi,Range::all()) += EigVecs(bi,bj)*Bands(bj,Range::all());
  
  Bands = RotBands;
  // And compute rotated HBands
  RotBands = 0.0;
  for (int bi=0; bi<NumBands; bi++)
    for (int bj=0; bj<NumBands; bj++)
      RotBands(bi,Range::all()) += EigVecs(bi,bj)*HBands(bj,Range::all());
  HBands = RotBands;

  zVec Tc(Bands.extent(1));
  for (int bi=0; bi<NumBands; bi++) {
    Tc = 0.0;
    c.reference (Bands(bi, Range::all()));
    H.Kinetic.Apply (c, Tc);
    CG.T(bi) = realconjdot (c, Tc);
  }
  // HACK HACK HACK
  // CG.ApplyH();

//   for (int bi=0; bi<NumBands; bi++) {
//     c.reference (Bands(bi,Range::all()));
//     Hc.reference (HBands(bi, Range::all()));
//     double E = realconjdot(c, Hc);
//     perr << "EigVal(bi) = " << EigVals(bi) << "cHc = " << E << endl;
//   }
}


void
MPISystemClass::SolveLDA()
{
  if (UseSubspaceRotation)
    CG.SetOrthoType (ORTHO_ALL);
  else
    CG.SetOrthoType (ORTHO_LOWER);

  // Make sure we have enough bands
  assert (NumBands >= (NumElecs+1)/2);
  int numSCIters = 0;
  bool SC = false;
  int highestOcc = (NumElecs+1)/2-1;
  double lastE = 500.0;
  ChargeMixer->Reset();
  
  double startLambda = ChargeMixer->GetLambda();

  while (!SC) {
    MixChargeDensity();
    CalcVHXC();
    CG.SetTolerance(1.0e-4);
    if (numSCIters == 0 && ConfigNum == 0) {
      CG.InitBands();
      CG.ApplyH();
      SubspaceRotate();
    }
    CG.ApplyH();
    if (UseSubspaceRotation) 
      SubspaceRotate();

    CG.Solve();
    SC = (fabs(CG.Energies(highestOcc)-lastE) < 1.0e-6);
    /// Now make sure we all k-points agree that we're
    /// self-consistent.
    if (BandComm.MyProc() == 0) 
      kComm.AllAnd(SC);
    BandComm.Broadcast(0, SC);

    lastE = CG.Energies(highestOcc);
    CalcOccupancies();
    CalcChargeDensity();
    numSCIters ++;
    /// Stability control:  if convergence is taking forever, we're
    /// probably  charge sloshing. Then increase the Kerker
    /// parameter
    if ((numSCIters % 50) == 0)
      ChargeMixer->SetLambda (1.2*ChargeMixer->GetLambda());
  }
  // Reset initial lambda
  ChargeMixer->SetLambda(startLambda);

//   if (BandComm.MyProc() == 0) {
//     if (kComm.MyProc() == 0) {
//       IOSectionClass out;
//       char fname[100];
//       snprintf (fname, 100, "VHXC%d.h5", ConfigNum);
//       out.NewFile (fname);
//       out.WriteVar("VH", VH);
//       out.WriteVar("VXC", VXC);
//       out.WriteVar("VHXC", VHXC);
//       out.WriteVar("Rho", Rho_r);
//       out.WriteVar ("E", CG.Energies);
//       LocalPotFFTClass &pot = *((LocalPotFFTClass*)H.Vion);
//       Array<double,3> magV(VH.shape());
//       for (int ix=0; ix<magV.extent(0); ix++)
// 	for (int iy=0; iy<magV.extent(1); iy++)
// 	  for (int iz=0; iz<magV.extent(2); iz++)
// 	    magV(ix,iy,iz) = mag(pot.Vr(ix,iy,iz));
//       out.WriteVar("Vion", magV);
//       out.CloseFile();
//     }
//   }
  Bands2 = Bands1;
  Rions2 = Rions1;
  Bands1 = Bands;
  Rions1 = Rions;
  Rho_r2 = Rho_r1;
  Rho_r1 = Rho_r;
  ConfigNum++;
}


double
Residue (Array<Vec3,1> &R0, Array<Vec3,1> &R1, Array<Vec3,1> &R2,
	 double alpha)
{
  double residue = 0;
  for (int i=0; i<R0.size(); i++)
    for (int dim=0; dim < 3; dim++) {
      double diff = R0(i)[dim] -(alpha*R1(i)[dim] + (1.0-alpha)*R2(i)[dim]);
      residue += diff*diff;
    }
  return residue;
}

/// See Dario Alfe, Computer Physics Communications 118, 31-33 (1999)
void
MPISystemClass::ExtrapolateCharge()
{


}

/// See Arias et al. Phys. Rev. B 45, 1538 (1992)
void 
MPISystemClass::DoMDExtrap()
{
  if (ConfigNum > 1) {
    double e0, e1, e2;
    e0  = Residue(Rions, Rions1, Rions2,  0.0);
    e1  = Residue(Rions, Rions1, Rions2,  1.0);
    e2  = Residue(Rions, Rions1, Rions2,  2.0);
    
    double a = 0.5*(e2-2.0*e1+e0);
    double b = 0.5*(4.0*e1 - e2 - 3.0*e0);
    double c = e0;
    double alpha;
    if (a > 0.0)  // we have a minimum somewhere
      alpha = -b/(2.0*a);
    else if (e0>e2)
      alpha = 2.0;
    else
      alpha = 0.0;
    
    alpha = max(0.0, alpha);
    alpha = min(2.0, alpha);

    cerr << "alpha = " << alpha << endl;
    Bands = alpha * Bands1 + (1.0-alpha)*Bands2;
    GramSchmidt(Bands);

    // Now extrapolate density
    Rho_r = alpha * Rho_r1 + (1.0-alpha)*Rho_r2;
    // For the benefit of the charge mixer
    NewRho = Rho_r;
  }
}


void
MPISystemClass::CalcOccupancies()
{
  ///////////////////////////////////////
  // Calculate mu from 0.5*(HOMO+LUMO) //
  ///////////////////////////////////////
  // First, gather all the energies from the k-Point clones.

  Array<double,2> allEnergies(Numk, NumBands), occ(Numk, NumBands);
  if (BandComm.MyProc() == 0) {
    allEnergies(Myk, Range::all()) = CG.Energies;
    kComm.AllGatherRows(allEnergies);
  }
  BandComm.Broadcast (0, allEnergies);
  // Now processor zero has all the energies;
  vector<double> Esort;
  for (int i=0; i<allEnergies.extent(0); i++)
    for (int j=0; j<allEnergies.extent(1); j++)
      Esort.push_back(allEnergies(i,j));
  sort (Esort.begin(), Esort.end());
  int HOMO = (NumElecs*Numk+1)/2-1;
  assert ((HOMO+1) < Esort.size());
  double muLo = Esort[0];
  double muHi = Esort[Esort.size()-1];
  double numkInv = 1.0/(double)Numk;
  double mu, occSum;
  while (fabs(occSum-(double)NumElecs) > 1.0e-10) {
    occSum = 0.0;
    mu = 0.5*(muHi + muLo);
    for (int ki=0; ki<Numk; ki++)
      for (int bi=0; bi<NumBands; bi++) {
	occ(ki, bi) = 2.0*numkInv*Smearer.S(allEnergies(ki, bi), mu);
	occSum += occ(ki, bi);
      }
    if (occSum < (double) NumElecs)
      muLo = mu;
    else if (occSum > (double) NumElecs)
      muHi = mu;
  }
  

//   int LUMO = (NumElecs*Numk)/2;
//   double mu = 0.5*(Esort[HOMO] + Esort[LUMO]);
//   /// Broadcast mu to all procs in kBand;
//   if (BandComm.MyProc() == 0)
//     kComm.Broadcast(0, mu);
//   /// Now proc 0 of each k-point clone has mu.
//   /// Now broadcast mu from each BandComm Proc 0 to the rests of the
//   /// band holders
//   BandComm.Broadcast (0, mu);
//   /// Now everybody has mu!  Let's do the Methfessel-Paxton occupation 
//   double totalOcc = 0.0;
//   for (int ki=0; ki<Numk; ki++)
//     for (int bi=0; bi<NumBands; bi++) {
//       occ(ki, bi) = Smearer.S(allEnergies(ki, bi), mu);
//       totalOcc += occ(ki, bi);
//     }
  /// Now normalize to to make sure we have exactly the right number
  /// of electrons 
  for (int bi=0; bi<NumBands; bi++) 
    Occupancies(bi) = (double)(NumElecs)/occSum * occ(Myk, bi);
  perr << "Occupancies = " << Occupancies << endl;
}

void
MPISystemClass::CalcChargeDensity()
{
  int nx, ny, nz;
  FFT.GetDims(nx,ny,nz);
  if (TempRho.shape() != TinyVector<int,3>(nx,ny,nz))
    TempRho.resize(nx,ny,nz);
  if (NewRho.shape() != TinyVector<int,3>(nx,ny,nz))
    NewRho.resize(nx,ny,nz);
  TempRho = 0.0;
  
  zVec band(Bands.extent(1));
  /// We assume that each of the bands is already normalized.
  for (int bi=CG.GetFirstBand(); bi<=CG.GetLastBand(); bi++) {
    band = Bands(bi, Range::all());
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
  double prefactor = volInv; // /(double)(nx*ny*nz);
  NewRho *= prefactor;
  double vCell = H.GVecs.GetBoxVol()/(double)(nx*ny*nz);
  double totalCharge = 0.0;
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	totalCharge += vCell * NewRho(ix,iy,iz);
  perr << "total NewRho charge = " << totalCharge << endl;
}



/// Calculate the hartree, exchange, and correlation potentials.
/// We assume that upon entry, Rho_r contains the charge density 
/// in real space and Rho_G contains it in reciporical space.
/// On exit, VH will contain the Hartree potential and VXC will
/// contain the exchange-correlation potential, both in real space.
void 
MPISystemClass::CalcVHXC()
{
  double totalCharge = 0.0;
  double vol = GVecs.GetBoxVol();
  double Ngrid = Rho_r.size();

  ///////////////////////
  // Hartree potential //
  ///////////////////////
  h_G = Rho_G;
  // Compute V_H in reciporical space
  double hartreeTerm = 0.0;
  //  double prefact = (4.0*M_PI/H.GVecs.GetBoxVol());
  double prefact = 4.0*M_PI;
  EH = 0.0;
  for (int i=0; i<h_G.size(); i++) {
    EH += norm(h_G(i)) * H.GVecs.DeltaGInv2(i);
    h_G(i) *= prefact*H.GVecs.DeltaGInv2(i);
  }
  EH *= 0.5*prefact*vol;
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
	EXC += exc*Rho_r(ix,iy,iz);
      }
  EXC *= vol/(double)Ngrid;
  perr<< "EH = " << EH << "   EXC = " << EXC << endl;
  VHXC = VH + VXC;
}


void
MPISystemClass::CalcRadialChargeDensity()
{
  perr << "Calculating atomic charge density:\n";
  double maxBox = max(Box[0], max(Box[1], Box[2]));
  double rmax = max (maxBox, 100.0);
  AtomGrid.Init (1.0, rmax);
  DFTAtom atom;
  atom.NewMix = 0.5;
  atom.RadialWFs.resize(1);
  atom.RadialWFs(0).n = 1;
  atom.RadialWFs(0).l = 0;
  atom.RadialWFs(0).Occupancy = 1.0;
  atom.RadialWFs(0).Energy = -0.15;
  atom.SetGrid (&AtomGrid);
  atom.SetBarePot (PH);
  atom.Solve();
  Array<double,1> rho(AtomGrid.NumPoints);
  for (int i=0; i<rho.size(); i++) {
    double r = AtomGrid(i);
    double u = atom.RadialWFs(0).u(i);
    rho(i) = u*u/(r*r);
  }
  RadialChargeDensity.Init (&AtomGrid, rho);
  perr << "Done atomic charge density.\n";
}

/// Initializes the density to a superposition of the atomic charge
/// densities from an Atomic DFT calculation.
void
MPISystemClass::InitNewRho()
{
  NewRho = 0.0;
  double nxInv = 1.0/(double)(Nx-1);
  double nyInv = 1.0/(double)(Ny-1);
  double nzInv = 1.0/(double)(Nz-1);
  Vec3 boxInv = Vec3(1.0/Box[0], 1.0/Box[1], 1.0/Box[2]);
  Vec3 r;
  double chargePerAtom = (double)NumElecs/(double)Rions.size();
  perr << "chargePerAtom = " << chargePerAtom << endl;
  double totalCharge = 0.0;
  double cellVol = Box[0]*Box[1]*Box[2]/(double)(Nx*Ny*Nz);
  for (int ix=0; ix<Nx; ix++) {
    r[0] = (nxInv*(double)ix-0.5)*Box[0];
    for (int iy=0; iy<Ny; iy++) {
      r[1] = (nyInv* (double)iy-0.5)*Box[1];
      for (int iz=0; iz<Nz; iz++) {
	r[2] = (nzInv*(double)iz-0.5)*Box[2];
	for (int ion=0; ion<Rions.size(); ion++) {
	  Vec3 disp = r-Rions(ion);
	  /// Find nearest image
	  disp[0] -= nearbyint(disp[0]*boxInv[0])*Box[0];
	  disp[1] -= nearbyint(disp[1]*boxInv[1])*Box[1];
	  disp[2] -= nearbyint(disp[2]*boxInv[2])*Box[2];
	  double dist = sqrt(dot(disp, disp));
	  double rho = RadialChargeDensity(dist);
	  if (rho < 0.0)
	    perr << "Negative rho at (ix,iy,iz) = (" 
		 << ix << ", " << iy << ", " << iz << ")\n";
	  NewRho(ix,iy,iz) += chargePerAtom*rho;
	}
	totalCharge += NewRho(ix,iy,iz);
      }
    }
  }
  totalCharge *= cellVol;
  NewRho *= (double)NumElecs/totalCharge;
}
	 
void 
MPISystemClass::MixChargeDensity()
{
  ChargeMixer->Mix(NewRho, Rho_r, Rho_G);
}



// double
// PlaneWavesMPI::CalcHartreeTerm(int band)
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


