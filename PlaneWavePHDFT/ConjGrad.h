#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include "Hamiltonian.h"

class ConjGrad
{
protected:
  Hamiltonian &H;
  zVec c;
  zVec cnext, Hc, Phi;
  zVec Phip, Phipp, Xi, Eta;
  complex<double> EtaXiLast;
  double E0;
  void CalcPhiSD();
  void CalcPhiCG();
  void Setup();
  bool IsSetup;
  void Precondition();
  double T;
  int iter, NumBands;
  int CurrentBand, LastBand;
  void PrintOverlaps();
public:
  Array<complex<double>,2> Bands;
  void Iterate(int band);
  ConjGrad (Hamiltonian &h, int numBands) : 
    H(h), IsSetup(false), EtaXiLast(0.0, 0.0), iter(0), NumBands(numBands),
    LastBand(-1)
  {
    // Do nothing for now
  }
};

#endif
