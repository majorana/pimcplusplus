#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "../FFT/FFT.h"

typedef TinyVector<int,3> Int3;
typedef Array<complex<double>,1> zVec;

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
  { return Index(i); }
  
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
  // Puts a linear c-vector into the fft box.
  FFTBox& operator= (const zVec &c);
  
  // Puts appropriate elements of fftbox into a c vector.
  void GetVec (zVec &c);
  
  FFTBox (GVecsClass &gvecs) : GVecs (gvecs)
  {
    int nx, ny, nz;
    gvecs.GetFFTBoxSize(nx,ny,nz);
    resize(nx, ny, nz);
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

class Hamiltonian : public HamiltonianBase
{
private:
  CoulombClass Coulomb;
  KineticClass Kinetic;
public:
  GVecsClass GVecs;
  void Apply (const zVec &c, zVec &Hc);
  Hamiltonian (Vec3 box, double kcut, double z) :
    Kinetic (GVecs), Coulomb(z, GVecs), HamiltonianBase(GVecs)
  {
    GVecs.Set (box, kcut);
  }
};


#endif
