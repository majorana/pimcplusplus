#include "Coulomb.h"

void
CoulombClass::Setup()
{

}


void 
CoulombClass::Vmatrix (Array<complex<double>,2> &vmat)
{
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<vmat.rows(); i++) 
    for (int j=0; j<=i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      complex<double> val = -4.0*volInv*s*M_PI*Z/dot(diff,diff);
      vmat(i,j) = val;
      vmat(j,i) = conj(val);
    }
  for (int i=0; i<vmat.rows(); i++)
    vmat(i,i) = 0.0;
}

void 
CoulombClass::Apply(const zVec &c, zVec &Vc)
{
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<c.size(); i++) {
    for (int j=0; j<i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      complex<double> val = -4.0*volInv*s*M_PI*Z/dot(diff,diff);
      Vc(i) += val*c(j);
      Vc(j) += conj(val)*c(i);
    }
  }
}


void
CoulombClass::SetIons (const Array<Vec3,1> &rions)
{
  if (Rions.size() != rions.size())
    Rions.resize(rions.size());
  Rions = rions;
}
