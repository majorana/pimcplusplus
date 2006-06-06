#ifndef PLANE_WAVES_MPI_H
#define PLANE_WAVES_MPI_H

#include "ConjGradMPI.h"
#include "../Splines/CubicSpline.h"
#include "../Atom/DFTAtom.h"

class CommunicatorClass;

class MPISystemClass
{
protected:
  Array<Vec3,1> Rions;
  FFTBox FFT;
  HamiltonianClass H;
  Array<complex<double>,2> Bands, LastBands;
  Array<double,1> Occupancies;
  ConjGradMPI CG;
  Vec3 Box;
  double kCut;
  int NumBands, NumElecs;
  Potential *PH;
  bool UseFFT;
  CommunicatorClass &BandComm, &kComm;
  bool MDExtrap, FirstTime;
  /////////////////////////
  /// LDA-related stuff ///
  /////////////////////////
  bool UseLDA;
  int Nx, Ny, Nz;
  // Array<complex<double>,3> Phip_r, Psi_r;
  void CalcOccupancies();
  void CalcChargeDensity();
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
  void InitRho_r();
  void SolveLDA();
  /// HACK HACK HACK
public:
  void InitLDA();
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
    : CG(H, Bands, bandcomm, kcomm, VHXC), 
      FFT(GVecs), H(GVecs, FFT), 
      NumBands(numBands), BandComm(bandcomm), kComm(kcomm), 
      MDExtrap(mdextrap), FirstTime(true), NumElecs(numElecs)
  {

  }
};

#endif
