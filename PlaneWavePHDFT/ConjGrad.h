#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include "Hamiltonian.h"

class ConjGrad
{
protected:
  Hamiltonian &H;
  zVec c, cnext, Hc, Hcnext, Phi, LastPhi;
  zVec Phip, Phipp, Xi, Xip, Eta, Etap;
  complex<double> EtaXiLast;
  double E0;
  void CalcPhiSD();
  void CalcPhiCG();
  void Setup();
  bool IsSetup;
  void Precondition();
  double T;
public:

  void Iterate();
  ConjGrad (Hamiltonian &h) : H(h), IsSetup(false),
			      EtaXiLast(0.0, 0.0)
  {
    // Do nothing for now
  }
};

#endif
