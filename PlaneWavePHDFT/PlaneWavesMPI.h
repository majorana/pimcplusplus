#ifndef PLANE_WAVES_MPI_H
#define PLANE_WAVES_MPI_H

#include "ConjGradMPI.h"
#include "../Splines/CubicSpline.h"
#include "../Atom/DFTAtom.h"
#include "FermiSmear.h"

class CommunicatorClass;

class MPISystemClass
{
protected:
  ///////////////////
  // Basic storage //
  ///////////////////
  /// The positions of the ions
  Array<Vec3,1> Rions;

  /// Simulation box size
  Vec3 Box;

  /// The reciporical-space cutoff
  double kCut;
  
  /// The number of bands and the number of electrons.  For technical
  /// reasons, it is often necessary to solve for unoccupied bands
  int NumBands, NumElecs;
  
  /// The electron-ion potential -- presently, this only works for
  /// systems containing one element type
  Potential *PH;

  /// Object for applying the Hamiltonian to a vector.
  HamiltonianClass H;

  // The band coefficients and the Hamiltonian applied to them
  Array<complex<double>,2> Bands, HBands;

  // Band occupancies
  Array<double,1> Occupancies;

  // The conjugate gradient solver
  ConjGradMPI CG;
  
  /// The workhorse for doing the FFTs.
  FFTBox FFT;

  /// The processors are divided into a number of groups, each of
  /// which work on a different k-point.  BandComm is an MPI
  /// communicator which includes all the processors working on this
  /// processor groups k-point.  kComm is a communicator which
  /// includes the root processor of each BandComm.  It is used to
  /// determine occupancies and sum the charge density over k-points.
  CommunicatorClass &BandComm, &kComm;
  /// The number of k-points and this processors k-point number
  int Numk, Myk;

  ///////////
  // Flags //
  ///////////
  bool UseFFT;
  bool UseMDExtrap, FirstTime;
  bool UseLDA;
  bool UseSubspaceRotation;

  ////////////////
  // IO related //
  ////////////////
  void Read (IOSectionClass &in);
  /// This stores the file location for storing any output
  IOSectionClass OutSection;

  //////////////////////
  // MD extrapolation //
  //////////////////////
  // The _1's contain the data from the previous ion configuration and
  // the _2's from the one before that.  We need both to extrapolate.
  Array<complex<double>,2> Bands1, Bands2;
  Array<double,3> Rho_r1, Rho_r2;
  Array<Vec3,1> Rions1, Rions2;

  ///////////////////////
  // Subspace rotation //
  ///////////////////////
  Array<complex<double>,2> Hmat, EigVecs, RotBands;
  Array<double,1> EigVals;
  void SubspaceRotate();
  /// Applies the current hamiltonian to all bands, storing in Hbands
  void ApplyH();

  /////////////////////////
  /// LDA-related stuff ///
  /////////////////////////
  int ConfigNum;
  int Nx, Ny, Nz;
  // Array<complex<double>,3> Phip_r, Psi_r;
  void CalcOccupancies();
  MethfesselPaxton Smearer;
  void CalcChargeDensity();
  ChargeMixerClass *ChargeMixer;
  void MixChargeDensity();
  /// These store the electron charge density in real space on the FFT
  /// grid.  TempRho is used to compute the charge density for my own
  /// k-point. 
  Array<double,3> TempRho, NewRho, Rho_r;
  zVec Rho_G;
  void CalcVHXC();
  /// The Hatree and exchange-correlation potentials in real space --
  /// they should always be real in real space
  Array<double,3> VH, VXC, VHXC;
  /// The components of the Hartree term in k-space
  zVec h_G;
  /// The Hartree and exchance-corelation energies
  double EH, EXC;
  /// This radial density is used to calculate 
  OptimalGrid AtomGrid;
  CubicSpline RadialChargeDensity;
  void CalcRadialChargeDensity();
  void InitNewRho();
  /// HACK HACK HACK
public:
  void InitLDA();
  void DoMDExtrap();
public:
  GVecsClass GVecs;
  void Setup(Vec3 box, Vec3 k, double kcut, Potential &ph, 
	     bool useLDA, bool useFFT=true);
  void Setup(Vec3 box, Vec3 k, double kcut, double z, 
	     bool useLDA, bool useFFT=true);
  void SetIons (const Array<Vec3,1> &rions);
  inline Vec3 GetIonPos(int i) { return Rions(i); }
  inline const Array<double,3>& GetDensity()
  { return Rho_r; }
  void Setk (Vec3 k);
  void DiagonalizeH();
  void SolveLDA();
  inline double GetEnergy(int band) { return CG.Energies(band); }
  inline int GetNumBands() { return NumBands; }

  /// Gets the FFT box dimensions.
  inline void GetBoxDims(int &nx, int &ny, int &nz)
  { FFT.GetDims(nx, ny, nz); }
  /// This FFT's the desired band into real space.
  void SetRealSpaceBandNum(int num);
  inline complex<double> RealSpaceBand (int ix, int iy, int iz)
  { return FFT.rBox(ix, iy, iz); }
  

  void CalcChargeDensity(Array<double,3> &rho);
  void WriteXSFFile(string filename);

  MPISystemClass(int numBands, int numElecs,
		 CommunicatorClass &bandcomm, CommunicatorClass &kcomm,
		 bool useLDA=false, bool mdextrap=false) 
    : CG(H, Bands, HBands, bandcomm, kcomm, VHXC), 
      FFT(GVecs), H(GVecs, FFT), 
      NumBands(numBands), BandComm(bandcomm), kComm(kcomm), 
      UseMDExtrap(mdextrap), FirstTime(true), NumElecs(numElecs),
      ConfigNum(0), UseSubspaceRotation(true)
  {

  }
};

#endif
