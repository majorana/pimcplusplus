#include "Hamiltonian.h"
#include "../PH/kSpacePH.h"

void GVecsClass::Set (Vec3 box, double kcut)
{
  Box = box;
  kBox[0]=2.0*M_PI/box[0]; kBox[1]=2.0*M_PI/box[1]; kBox[2]=2.0*M_PI/box[2];

  int maxX, maxY, maxZ;
  maxX = (int) ceil(kcut/kBox[0]);
  maxY = (int) ceil(kcut/kBox[1]);
  maxZ = (int) ceil(kcut/kBox[2]);

  // The FFT box must be twice the size as the maximum G in each direction.
  Nx = 2*(2*maxX+1);
  Ny = 2*(2*maxY+1);
  Nz = 2*(2*maxZ+1);

  int numG = 0;
  Vec3 G;
  /// First, count k-vectors
  for (int ix=-maxX; ix<=maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-maxY; iy<=maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (kcut*kcut))
	  numG++;
      }
    }
  }
  GVecs.resize(numG);
  Indices.resize(numG);
  numG = 0;
  /// Now assign
  for (int ix=-maxX; ix<=maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-maxY; iy<=maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (kcut*kcut)) {
	  GVecs(numG) = G;
	  Indices(numG) = Int3(ix,iy,iz);
	  numG++;
	}
      }
    }
  }
  cerr << "Using " << numG << " G-vectors.\n";
}


FFTBox &
FFTBox::operator=(const Array<complex<double>,1> &c)
{


}


void KineticClass::Setup()
{
  halfG2.resize(GVecs.size());
  for (int i=0; i<GVecs.size(); i++)
    halfG2(i) = 0.5*dot(GVecs(i),GVecs(i));
  IsSetup = true;
}

void KineticClass::Apply(const zVec &c, zVec &Kc)
{
  if (!IsSetup)
    Setup();
  for (int i=0; i<c.size(); i++)
    Kc(i) += halfG2(i)*c(i);
}


void CoulombClass::Apply(const zVec &c, zVec &Vc)
{
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<c.size(); i++) 
    for (int j=0; j<c.size(); j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      Vc(i) += -4.0*volInv*M_PI*Z/dot(diff,diff)*c(j);
    }
}


void Hamiltonian::Apply(const zVec &c, zVec &Hc)
{
  Hc = 0.0;
  Kinetic.Apply (c, Hc);
  Coulomb.Apply (c, Hc);
}
