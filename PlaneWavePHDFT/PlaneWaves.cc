#include "PlaneWaves.h"

void
SystemClass::Setup (Vec3 box, Vec3 k, double kcut, Potential &ph,
		    bool useFFT)
{
  Box = box;
  GVecs.Set (box, k, kcut);
  H.SetIonPot (ph, useFFT);
  Bands.resize (NumBands, GVecs.size());
}


void 
SystemClass::Setup (Vec3 box, Vec3 k, double kcut, double z,
		    bool useFFT)
{
  Box = box;
  GVecs.Set (box, k, kcut);
  H.SetIonPot (z, useFFT);
  H.Setk(k);
  cerr << "NumBands = " << NumBands << endl;
  Bands.resize (NumBands, GVecs.size());
}

void
SystemClass::SetIons (const Array<Vec3,1> &rions)
{
  if (Rions.size() != rions.size())
    Rions.resize(rions.size());
  Rions = rions;
  H.SetIons(rions);
}


void 
SystemClass::DiagonalizeH ()
{
  for (int i=0; i<Bands.rows(); i++) {
    CG.Solve(i);
    fprintf (stderr, "Energy(%d) = %15.12f\n", i, CG.Energies(i)* 27.211383);
  }
}
