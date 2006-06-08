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
  Array<complex<double>,2> &HBands;
  Array<complex<double>,2> lastPhis;
  double Tolerance;
  CommunicatorClass &BandComm, &kComm;
  /// Returns the maximum residual of all the bands
  double Iterate();
  void CollectBands();
  int MyFirstBand, MyLastBand;
  Array<double,1> Residuals;
  void CheckOverlaps();
  /// For LDA
  Array<double,3> &VHXC;

public:
  void Setup();

  /// Stores the band eigenenergies
  Array<double,1> Energies;
  void InitBands();
  void Solve();
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  inline int GetFirstBand() { return MyFirstBand; }
  inline int GetLastBand () { return MyLastBand; }

  ConjGradMPI (HamiltonianClass &h,
	       Array<complex<double>,2> &bands,
	       Array<complex<double>,2> &hbands,
	       CommunicatorClass &bandComm, 
	       CommunicatorClass &kcomm,
	       Array<double,3> &vhxc) : 
    H(h), IsSetup(false), iter(0),
    LastBand(-1), Bands(bands), HBands(hbands), Tolerance (1.0e-6),
    BandComm(bandComm), kComm(kcomm), VHXC(vhxc)
  {
    // Do nothing for now
  }
};

#endif
