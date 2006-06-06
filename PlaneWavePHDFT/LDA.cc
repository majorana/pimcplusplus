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
  TempRho.resize (Nx, Ny, Nz);
  VH.resize      (Nx, Ny, Nz);
  VXC.resize     (Nx, Ny, Nz);
  VHXC.resize    (Nx, Ny, Nz);
  Rho_r.resize   (Nx, Ny, Nz);
  h_G.resize     (NDelta);
  Rho_G.resize   (NDelta);
  CalcRadialChargeDensity();
  Smearer.SetOrder(2);
  Smearer.SetWidth(0.05);
}


void
MPISystemClass::SolveLDA()
{
  int numSCIters = 0;
  bool SC = false;
  while (!SC) {
    MixChargeDensity();
    CalcVHXC();
    IOSectionClass out;
    out.NewFile ("VHXC.h5");
    out.WriteVar("VH", VH);
    out.WriteVar("VXC", VXC);
    out.WriteVar("VHXC", VHXC);
    out.CloseFile();
    cerr << "After CaclVHXC.\n";
    CG.SetTolerance(1.0e-6);
    cerr << "Before CG.Solve.\n";
    CG.Solve();
    cerr << "Before CalcOccupancies.\n";
    CalcOccupancies();
    cerr << "Before CalcChargeDensity.\n";
    CalcChargeDensity();
    numSCIters ++;
  }

}

void
MPISystemClass::CalcOccupancies()
{
  ///////////////////////////////////////
  // Calculate mu from 0.5*(HOMO+LUMO) //
  ///////////////////////////////////////
  // First, gather all the energies from the k-Point clones.
  int numK = kComm.NumProcs();
  int myK  = kComm.MyProc();
  Array<double,2> allEnergies(numK, NumBands), occ(numK, NumBands);
  if (BandComm.MyProc() == 0) {
    allEnergies(myK, Range::all()) = CG.Energies;
    kComm.AllGatherRows(allEnergies);
  }
  BandComm.Broadcast (0, allEnergies);
  // Now processor zero has all the energies;
  vector<double> Esort;
  for (int i=0; i<allEnergies.extent(0); i++)
    for (int j=0; j<allEnergies.extent(1); j++)
      Esort.push_back(allEnergies(i,j));
  sort (Esort.begin(), Esort.end());
  int HOMO = (NumElecs*numK+1)/2-1;
  int LUMO = (NumElecs*numK)/2;
  double mu = 0.5*(Esort[HOMO] + Esort[LUMO]);
  /// Broadcast mu to all procs in kBand;
  if (BandComm.MyProc() == 0)
    kComm.Broadcast(0, mu);
  /// Now proc 0 of each k-point clone has mu.
  /// Now broadcast mu from each BandComm Proc 0 to the rests of the
  /// band holders
  BandComm.Broadcast (0, mu);
  /// Now everybody has mu!  Let's do the Methfessel-Paxton occupation 
  double totalOcc;
  for (int ki=0; ki<numK; ki++)
    for (int bi=0; bi<Bands.size(); bi++) {
      occ(ki, bi) = Smearer.S(CG.Energies(bi), mu);
      totalOcc += occ(ki, bi);
    }
  /// Now normalize to to make sure we have exactly the right number
  /// of electrons 
  for (int bi=0; bi<Bands.size(); bi++) 
    Occupancies(bi) = (double)(NumElecs*numK)/totalOcc * occ(myK, bi);
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
  zVec band;

  /// We assume that each of the bands is already normalized.
  for (int bi=CG.GetFirstBand(); bi<CG.GetLastBand(); bi++) {
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
  double prefactor = volInv/(double)(nx*ny*nz);
  NewRho *= prefactor;
  double vCell = H.GVecs.GetBoxVol()/(double)(nx*ny*nz);
  double totalCharge = 0.0;
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	totalCharge += vCell * NewRho(ix,iy,iz);
  cerr << "totalCharge = " << totalCharge << endl;
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
    EH += norm(h_G(i)) * H.GVecs.DeltaGInv2(i);
    h_G(i) *= prefact*H.GVecs.DeltaGInv2(i);
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
  VHXC = VH + VXC;
}


void
MPISystemClass::CalcRadialChargeDensity()
{
  cerr << "Calculating atomic charge density:\n";
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
  //atom.SolveInit();
  atom.Solve();
  Array<double,1> rho(AtomGrid.NumPoints);
  for (int i=0; i<rho.size(); i++) {
    double r = AtomGrid(i);
    double u = atom.RadialWFs(0).u(i);
    rho(i) = u*u/(r*r);
  }
  RadialChargeDensity.Init (&AtomGrid, rho);
  cerr << "Done atomic charge density.\n";
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
	  NewRho(ix,iy,iz) += chargePerAtom*RadialChargeDensity(dist);
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


