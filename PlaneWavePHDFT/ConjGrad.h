#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include "Hamiltonian.h"

class ConjGrad
{
protected:
  Hamiltonian &H;
  zVec c, cnext, Hc, Hcnext, SDvec;
  double E0;
  void CalcSD();
  void Setup();
  bool IsSetup;
public:

  void Iterate();
  ConjGrad (Hamiltonian &h) : H(h), IsSetup(false)
  {
    // Do nothing for now
  }
};

#endif
