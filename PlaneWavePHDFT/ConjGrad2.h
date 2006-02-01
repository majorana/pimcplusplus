#ifndef CONJ_GRAD2_H
#define CONJ_GRAD2_H

#include "Hamiltonians.h"

class ConjGrad
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

public:
  void Setup();
  Array<double,1> Energies;
  void InitBands();
  void Solve(int band);
  inline void SetTolerance(double tol) { Tolerance = tol;}
  void PrintOverlaps();

  ConjGrad (HamiltonianClass &h, Array<complex<double>,2> &bands) : 
    H(h), IsSetup(false), EtaXiLast(0.0, 0.0), iter(0),
    LastBand(-1), Bands(bands), Tolerance (1.0e-8)
  {
    // Do nothing for now
  }
};

#endif
