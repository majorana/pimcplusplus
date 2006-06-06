#include "Hamiltonians.h"

void
HamiltonianClass::SetIonPot (double z, bool useFFT)
{
  if (useFFT)
    Vion = new CoulombFFTClass(z, GVecs, FFT);
  else
    Vion = new CoulombClass (z, GVecs);
}


void 
HamiltonianClass::SetIonPot (Potential &pot, bool useFFT)
{
  if (useFFT) {
    if (pot.IsPH())
      Vion = new PHPotFFTClass (pot, GVecs, FFT);
    else
      Vion = new LocalPotFFTClass (pot, GVecs, FFT);
  }
  else
    Vion = new PHPotClass (pot, GVecs);
}


void 
HamiltonianClass::SetIons (const Array<Vec3,1>& rions)
{
  Vion->SetIons(rions);
}

void
HamiltonianClass::Setk (Vec3 k)
{
  kPoint = k;
  Kinetic.Setk(k);
  Vion->Setk(k);
}

void
HamiltonianClass::Apply (const zVec &c, zVec &Hc)
{
  Hc = 0.0;
  Kinetic.Apply(c, Hc);
  Vion->Apply (c, Hc);
}

void
HamiltonianClass::Apply (const zVec &c, zVec &Hc,
			 Array<double,3> &VHXC)
{
  Hc = 0.0;
  Kinetic.Apply(c, Hc);
  Vion->Apply (c, Hc, VHXC);
}
