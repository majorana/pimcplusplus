#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "../FFT/FFT.h"
#include "VectorOps.h"


class GVecsClass
{
protected:
  Array<Vec3,1> GVecs;
  Array<Int3,1> Indices;
  Vec3 Box, kBox;
  int Nx, Ny, Nz;
  double kCut;

public:
  void Set (Vec3 box, double kcut);

  inline Vec3 operator()(int i) const
  { return GVecs(i); }

  inline Int3 Index(int i) const
  { return Indices(i); }
  
  inline int size()
  { return GVecs.size(); }

  inline double GetBoxVol() 
  { return Box[0]*Box[1]*Box[2]; }

  void GetFFTBoxSize (int &nx, int &ny, int &nz);
  void PutVecInFFTBox (zVec &c, 
		       FFT3D &fft);
  void GetVecFromFFTBox (zVec &c,
			 FFT3D &fft);
  
};

class FFTBox : public FFT3D
{
protected:
  GVecsClass &GVecs;
public:
  int Nx, Ny, Nz;

  // Puts a linear c-vector into the fft box.
  void PutkVec (const zVec &c);
  
  // Puts appropriate elements of fftbox into a c vector.
  void GetkVec (zVec &c);
  
  // Adds the appropriate elements of the fftbox to a c vector.
  void AddFromVec (const zVec &c);

  // Add the c vector to the contents of the fftbox.
  void AddToVec (zVec &c);
  
  void Setup()
  {
    GVecs.GetFFTBoxSize(Nx, Ny, Nz);
    resize(Nx, Ny, Nz);
  }

  FFTBox (GVecsClass &gvecs) : GVecs (gvecs)
  {
  }
};


class HamiltonianBase
{
protected:
  GVecsClass &GVecs;
public:
  // Adds H*c to the Hc vector.  Does not zero Hc before accumulating
  virtual void Apply (const zVec &c, zVec &Hc) = 0;
  HamiltonianBase (GVecsClass &gvecs) : GVecs(gvecs)
  {
    // do nothing for now
  }
};


class KineticClass : public HamiltonianBase
{
private:
  Array<double,1> halfG2;
  void Setup();
  bool IsSetup;
public:
  void Apply (const zVec &c, zVec &Kc);
  KineticClass (GVecsClass &gvecs) : 
    HamiltonianBase (gvecs),  IsSetup(false)
  {
  }
};


class CoulombClass : public HamiltonianBase
{
private:
  double Z;
public:
  void Apply (const zVec &c, zVec &Hc);
  CoulombClass (double z, GVecsClass &gvecs) : 
    HamiltonianBase (gvecs), Z(z)
  {
    // nothing for now
  }
};

class CoulombFFTClass : public HamiltonianBase
{
private:
  double Z;
  void Setup();
  bool IsSetup;
  Array<complex<double>,3> Vr;
public:
  FFTBox FFT;
  void Apply (const zVec &c, zVec &Hc);
  CoulombFFTClass (double z, GVecsClass &gvecs) : 
    HamiltonianBase (gvecs), Z(z), FFT(gvecs), IsSetup(false)
  {
    // nothing for now
  }
};

class Hamiltonian : public HamiltonianBase
{
private:
  CoulombClass Coulomb;
  CoulombFFTClass CoulombFFT;
  KineticClass Kinetic;
public:
  GVecsClass GVecs;
  void Apply (const zVec &c, zVec &Hc);
  Hamiltonian (Vec3 box, double kcut, double z) :
    Kinetic (GVecs), Coulomb(z, GVecs), HamiltonianBase(GVecs),
    CoulombFFT(z, GVecs)
  {
    GVecs.Set (box, kcut);
  }
};


#endif
