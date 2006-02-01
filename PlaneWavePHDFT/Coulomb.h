#ifndef COULOMB_H
#define COULOMB_H

#include "HamiltonianBase.h"

class CoulombClass : public VionBase
{
private:
  double Z;
public:
  void Setup();
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);

  CoulombClass (double z, GVecsClass &gvecs) : 
    VionBase (gvecs), Z(z)
  {
    // nothing for now
  }
};

#endif
