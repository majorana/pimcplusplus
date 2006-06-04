#include "HamiltonianBase.h"

void 
PWKineticClass::Setup()
{
  halfG2.resize(GVecs.size());
  for (int i=0; i<GVecs.size(); i++) {
    Vec3 Gpk = GVecs(i)+kPoint;
    halfG2(i) = 0.5*dot(Gpk, Gpk);
  }
  IsSetup = true;
}

void 
PWKineticClass::Setk (Vec3 k)
{
  kPoint = k;
  Setup();
}

void 
PWKineticClass::Apply(const zVec &c, zVec &Kc)
{
  if (!IsSetup)
    Setup();
  for (int i=0; i<c.size(); i++)
    Kc(i) += halfG2(i)*c(i);
}


void
VionBase::SetIons(const Array<Vec3,1> &rions)
{
  Rions.resize(rions.size());
  Rions = rions;
  if (StructureFactor.size() != GVecs.DeltaSize())
    StructureFactor.resize(GVecs.DeltaSize());
  StructureFactor = 0.0;
  for (int gi=0; gi<GVecs.DeltaSize(); gi++) {
    double c, s, phase;
    for (int i=0; i<rions.size(); i++) {
      // HACK HACK HACK -- trying minus sign
      phase = dot (GVecs.DeltaG(gi), rions(i));
      sincos (phase, &s, &c);
      StructureFactor(gi) += complex<double>(c, s);
    }
  }
}


void
VionBase::Apply   (const zVec &c, zVec &Hc, 
		   Array<complex<double>,3> &VHXC)
{
  cerr << "Applying VHXC potential not supported by non-FFT versions"
       << " of the potentials.\n"; 
  abort();
}
