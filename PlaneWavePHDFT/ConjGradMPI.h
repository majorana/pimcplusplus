#ifndef CONJ_GRAD_MPI_H
#define CONJ_GRAD_MPI_H

#include "Hamiltonian2.h"

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
  Array<complex<double>,2> SD, lastPhis, Phips, Hcs;
  double Tolerance;
  CommunicatorClass &Communicator;
  /// Returns the maximum residual of all the bands
  double Iterate();
  void CollectBands();
  void CalcSDs();
  int MyFirstBand, MyLastBand;
  Array<double,1> Residuals;
  void CheckOverlaps();
public:
  void Setup();
  Array<double,1> Energies;
  void InitBands();
  void Solve();
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  ConjGradMPI (HamiltonianClass &h, Array<complex<double>,2> &bands,
	       CommunicatorClass &comm) : 
    H(h), IsSetup(false), iter(0),
    LastBand(-1), Bands(bands), Tolerance (1.0e-8),
    Communicator(comm)
  {
    // Do nothing for now
  }
};

#endif
