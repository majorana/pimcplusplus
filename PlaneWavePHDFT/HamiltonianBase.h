#ifndef HAMILTONIAN_BASE_H
#define HAMILTONIAN_BASE_H

#include "VectorOps.h"
#include "../PH/kSpacePH.h"
#include "FFTBox.h"


class PWKineticClass
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

  PWKineticClass (GVecsClass &gvecs) : 
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
  virtual void Vmatrix (Array<complex<double>,2> &vmat) = 0;
  
  virtual void Setk    (Vec3 k) 
  { kPoint = k; }
  
  VionBase (GVecsClass &gvecs) 
    : GVecs(gvecs), IsSetup(false), kPoint (0.0, 0.0, 0.0)
  {
    // do nothing for now
  }
};


class HamiltonianClass
{
private:
public:
  FFTBox &FFT;
  GVecsClass &GVecs;
  PWKineticClass Kinetic;
  VionBase *Vion;
  Vec3 kPoint;

  void Apply (const zVec &c, zVec &Hc);
  
  void SetIonPot (double z, bool useFFT);
  void SetIonPot (Potential &ph, bool useFFT);
  void SetIons (const Array<Vec3,1>& rions);
  void Setk (Vec3 k);

  HamiltonianClass (GVecsClass &gvecs, FFTBox &fft) 
    : GVecs(gvecs), Kinetic(gvecs), Vion(NULL), FFT(fft)
  {

  }
};



#endif
