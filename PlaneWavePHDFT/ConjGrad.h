#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include "Hamiltonian.h"

inline void Normalize (zVec &c)
{
  double norm = 0.0;
  for (int i=0; i<c.size(); i++)
    norm += c(i).real()*c(i).real() + c(i).imag()*c(i).imag();
  norm = 1.0/sqrt(norm);
  for (int i=0; i<c.size(); i++)
    c(i) *= norm;
}

inline complex<double> conjdot(zVec &cA, zVec &cB)
{
  complex<double> z(0.0, 0.0);
  for (int i; i<cA.size(); i++)
    z += conj(cA(i))*cB(i);
  return z;
}

inline double realconjdot(zVec &cA, zVec &cB)
{
  double re = 0.0;
  for (int i; i<cA.size(); i++)
    re += real(conj(cA(i))*cB(i));
  return re;
}

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
