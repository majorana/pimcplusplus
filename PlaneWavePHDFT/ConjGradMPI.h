#ifndef CONJ_GRAD_MPI_H
#define CONJ_GRAD_MPI_H

#include "Hamiltonians.h"

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
  double Tolerance;
  CommunicatorClass &BandComm, &kComm;
  /// Returns the maximum residual of all the bands
  double Iterate();
  void CollectBands();
  int MyFirstBand, MyLastBand;
  Array<double,1> Residuals;
  void CheckOverlaps();
  /////////////////////////
  /// DFT-related stuff ///
  /////////////////////////
  FFTBox &FFT;
  bool UseLDA;
  /// Stores the electron charge density in real space on the FFT
  /// grid. 
  int Nx, Ny, Nz;
  Array<complex<double>,3> Phip_r, Psi_r;
  /// The components of the Hatree term in k-space
  zVec h_G;
  void CalcChargeDensity();
  /// MyRho is used to compute the charge density for my own k-point.  
  Array<double,3> TempRho, NewRho, Rho;
  void CalcVHXC();
  double CalcHartreeTerm(int band);
  double CalcXCTerm(int band);
  /// The Hatree and exchange-correlation potentials in real space
  Array<complex<double>,3> VH, VXC;
  /// The Hartree and exchance-corelation energies
  double EH, EXC;
public:
  void Setup();
  /// Stores the band eigenenergies
  Array<double,1> Energies;
  void InitBands();
  void Solve();
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  ConjGradMPI (HamiltonianClass &h, Array<complex<double>,2> &bands,
	       CommunicatorClass &bandComm, 
	       CommunicatorClass &kcomm,
	       FFTBox &fft) : 
    H(h), IsSetup(false), iter(0),
    LastBand(-1), Bands(bands), Tolerance (1.0e-6),
    BandComm(bandComm), kComm(kcomm), FFT(fft), UseLDA(false)
  {
    // Do nothing for now
  }
};

#endif
