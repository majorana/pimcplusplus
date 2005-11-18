#include "Forces.h"
#include "../PathDataClass.h"

void
ForcesClass::SetSpecies(int speciesNum)
{
  SpeciesNum = speciesNum;
  Fsum.resize(PathData.Path.Species(speciesNum).NumParticles);
  Fsum = 0.0;
  Counts = 0;
}

void
ForcesClass::Accumulate()
{
  PathClass &Path = PathData.Path;
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;  

  //////////////////////
  // Short-range part //
  //////////////////////
  for (int pi=0; pi<species.NumParticles; pi++) {
    int ptcl = pi+first;
    for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      PairActionFitClass &PA = 
	*(PathData.Actions.PairMatrix(SpeciesNum, species2));
      if (ptcl != ptcl2) 
	for (int link=0; link<(Path.NumTimeSlices()-1); link++) {
	  dVec r, rp;
	  double rmag, rpmag, du_dq, du_dz;
	  Path.DistDisp(link,link+1,ptcl, ptcl2, rmag, rpmag, r, rp);
	  double q = 0.5*(rmag+rpmag);
	  double z = 0.5*(rmag-rpmag);
	  double s2 = dot (r-rp, r-rp);
	  PA.Derivs(q,z,s2,0,du_dq, du_dz);
	  Vec3 rhat  = (1.0/rmag)*r;
	  Vec3 rphat = (1.0/rpmag)*rp;
	  Fsum(pi) += 0.5*(du_dq*(rhat+rphat) + du_dz*(rhat-rphat));
	  /// Now, subtract off long-range part that shouldn't be in
	  /// here 
	  if (PA.IsLongRange())
	    Fsum(pi) -= 0.5*(PA.Ulong(0).Deriv(rmag)*rhat+
			     PA.Ulong(0).Deriv(rpmag)*rphat);
	}
    }
  }

      
  //////////////////////
  // Long-range part  //
  //////////////////////





  Counts++;
}
