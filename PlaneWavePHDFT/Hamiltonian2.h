#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "VectorOps.h"
#include "../PH/kSpacePH.h"
#include "FFTBox.h"


class KineticClass
{
private:
  GVecsClass &GVecs;
  Array<double,1> halfG2;
  void Setup();
  Vec3 kPoint;
  bool IsSetup;
public:
  void Apply (const zVec &c, zVec &Kc);
  void Setk  (Vec3 k);

  KineticClass (GVecsClass &gvecs) : 
    GVecs (gvecs), IsSetup(false), kPoint(0.0, 0.0, 0.0)
  {
  }
};



class VionBase
{
protected:
  GVecsClass &GVecs;
  bool IsSetup;
  Vec3 kPoint;
  zVec StructureFactor;
  Array<Vec3,1> Rions;
public:
  // Adds H*c to the Hc vector.  Does not zero Hc before accumulating
  virtual void Apply   (const zVec &c, zVec &Hc) = 0;
  virtual void SetIons (const Array<Vec3,1> &rions); 
  virtual void Setup() = 0;
  
  virtual void Setk    (Vec3 k) 
  { kPoint = k; }
  
  VionBase (GVecsClass &gvecs) 
    : GVecs(gvecs), IsSetup(false), kPoint (0.0, 0.0, 0.0)
  {
    // do nothing for now
  }
};


class CoulombClass : public VionBase
{
private:
  double Z;
public:
  void Setup();
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);

  CoulombClass (double z, GVecsClass &gvecs) : 
    VionBase (gvecs), Z(z)
  {
    // nothing for now
  }
};

class CoulombFFTClass : public VionBase
{
private:
  double Z;
  Array<complex<double>,3> Vr;

  void Setup();
  void SetVr();
  FFTBox FFT;
public:
  void Apply (const zVec &c, zVec &Hc);
  void SetIons (const Array<Vec3,1> &rions);

  CoulombFFTClass (double z, GVecsClass &gvecs) : 
    VionBase (gvecs), Z(z), FFT(gvecs)
  {
    // nothing for now
  }
};

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
  void Setk(Vec3 k);
  
  PHPotClass (Potential &ph, GVecsClass &gvecs) :
    VionBase (gvecs), kPH(ph), VmatIsSet(false),
    SFIsSet(false)
  {

  }
};


class PHPotFFTClass : public VionBase
{
private:
  kSpacePH kPH;
  Array<Mat3,1> Fk;
  Array<cMat3,3> Fr;
  Array<double,1> Vk;
  Array<complex<double>,3> Vr;
  zVec Vc;
  zVecVec Gc;
  FFTBox      cFFT;
  FFTVecBox VecFFT;
  FFTMatBox MatFFT;
  void SetupkPotentials();
  void SetuprPotentials();
public:
  void Setup();
  void SetIons (const Array<Vec3, 1>& rions);

  void Apply (const zVec &c, zVec &Hc);
  PHPotFFTClass (Potential &ph, GVecsClass &gvecs) :
    VionBase (gvecs), kPH(ph),
    cFFT(gvecs), VecFFT(gvecs), MatFFT(gvecs)
  {

  }
};


class HamiltonianClass
{
private:
public:
  GVecsClass &GVecs;
  KineticClass Kinetic;
  VionBase *Vion;

  void Apply (const zVec &c, zVec &Hc);
  
  void SetIonPot (double z, bool useFFT);
  void SetIonPot (Potential &ph, bool useFFT);
  void SetIons (const Array<Vec3,1>& rions);

  HamiltonianClass (GVecsClass &gvecs) 
    : GVecs(gvecs), Kinetic(gvecs), Vion(NULL)
  {

  }
};


#endif
