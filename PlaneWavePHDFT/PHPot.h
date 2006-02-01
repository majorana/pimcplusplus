#ifndef PH_POT_H
#define PH_POT_H

#include "HamiltonianBase.h"

class PHPotClass : public VionBase
{
private:
  kSpacePH kPH;
  void Setup();
  Array<complex<double>,2> VGGp;
  Array<complex<double>,2> StructFact;
  void CalcStructFact();
  bool VmatIsSet, SFIsSet;
  void SetVmat();
public:
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);
  
  PHPotClass (Potential &ph, GVecsClass &gvecs) :
    VionBase (gvecs), kPH(ph), VmatIsSet(false),
    SFIsSet(false)
  {

  }
};


#endif
