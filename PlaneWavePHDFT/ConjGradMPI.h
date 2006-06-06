#ifndef CONJ_GRAD_MPI_H
#define CONJ_GRAD_MPI_H

#include "Hamiltonians.h"
#include "ChargeMixer.h"

class CommunicatorClass;

class ConjGradMPI
{
protected:
  HamiltonianClass &H;
  zVec c;
  zVec cnext, Hc, Phi;
  zVec Phip, Phipp, Xi, Eta;
  Array<complex<double>,1> EtaXiLast;
  double E0;
  double CalcPhiSD();
  double CalcPhiCG();
  bool IsSetup;
  void Precondition();
  Array<double,1> T;
  int iter;
  int CurrentBand, LastBand;
  Array<complex<double>,2> &Bands;
  Array<complex<double>,2> lastPhis;
  Array<double,1> Occupancies;
  double Tolerance;
  CommunicatorClass &BandComm, &kComm;
  /// Returns the maximum residual of all the bands
  double Iterate();
  void CollectBands();
  int MyFirstBand, MyLastBand;
  Array<double,1> Residuals;
  void CheckOverlaps();
  /////////////////////////
  /// LDA-related stuff ///
  /////////////////////////
  FFTBox &FFT;
  bool UseLDA;
  /// Stores the electron charge density in real space on the FFT
  /// grid. 
  int Nx, Ny, Nz;
  // Array<complex<double>,3> Phip_r, Psi_r;
  /// The components of the Hartree term in k-space
  zVec h_G;
  void CalcOccupancies();
  void CalcChargeDensity();
  void MixChargeDensity();
  /// MyRho is used to compute the charge density for my own k-point.  
  Array<double,3> TempRho, NewRho, Rho_r;
  zVec Rho_G;
  void CalcVHXC();
  /// The Hatree and exchange-correlation potentials in real space --
  /// they should always be real in real space
  Array<double,3> VH, VXC, VHXC;
  /// The Hartree and exchance-corelation energies
  double EH, EXC;
  int NumElecs;
public:
  void Setup();
  /// This sets the charge density.  Rho should have the same
  /// dimensions as the FFTBox.
  void SetChargeDensity(Array<double,3> &rho);
  /// This returns the current charge density
  void GetChargeDensity(Array<double,3> &rho);

  /// Stores the band eigenenergies
  Array<double,1> Energies;
  void InitBands();
  void Solve();
  void SolveLDA();
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  ConjGradMPI (HamiltonianClass &h, int numElecs,
	       Array<complex<double>,2> &bands,
	       CommunicatorClass &bandComm, 
	       CommunicatorClass &kcomm,
	       FFTBox &fft, bool useLDA=false) : 
    H(h), NumElecs(numElecs), IsSetup(false), iter(0),
    LastBand(-1), Bands(bands), Tolerance (1.0e-6),
    BandComm(bandComm), kComm(kcomm), FFT(fft), UseLDA(useLDA)
  {
    // Do nothing for now
  }
};

#endif
