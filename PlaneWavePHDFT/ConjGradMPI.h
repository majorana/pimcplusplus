#ifndef CONJ_GRAD2_H
#define CONJ_GRAD2_H

#include "Hamiltonian2.h"

class CommunicatorClass;

class ConjGradMPI
{
protected:
  HamiltonianClass &H;
  zVec c;
  zVec cnext, Hc, Phi;
  zVec Phip, Phipp, Xi, Eta;
  complex<double> EtaXiLast;
  double E0;
  double CalcPhiSD();
  double CalcPhiCG();
  bool IsSetup;
  void Precondition();
  double T;
  int iter;
  int CurrentBand, LastBand;
  Array<complex<double>,2> &Bands;
  double Tolerance;
  CommunicatorClass &Communicator;
  void CollectBands();
  int MyFirstBand, MyLastBand;
public:
  void Setup();
  Array<double,1> Energies;
  void InitBands();
  void Solve();
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  ConjGradMPI (HamiltonianClass &h, Array<complex<double>,2> &bands,
	    CommunicatorClass &comm) : 
    H(h), IsSetup(false), EtaXiLast(0.0, 0.0), iter(0),
    LastBand(-1), Bands(bands), Tolerance (1.0e-8),
    Communicator(comm)
  {
    // Do nothing for now
  }
};

#endif
