#include "ShortRangeClass.h"

double ShortRangeClass::Evaluate (int slice1, int slice2,
				  const Array<double,1> &activeParticles
				  int level)
{
  // First, sum the pair actions
  for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
    Path.DoPtcl(ptcl)=true;

  double TotalU = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = tau* (1<<level);
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = changedParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++) {
      if (Path.DoPtcl(ptcl2)){
	int PairIndex = PairMatrix(species1,
				   Path.ParticleSpeciesNum(ptcl2));

	for (int slice=startSlice;slice<endSlice;slice+=skip){
	  dVec r, rp;
	  double rmag, rpmag;

	  PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				 rmag, rpmag, r, rp);

	  double s2 = dot (r-rp, r-rp);
	  double q = 0.5 * (rmag + rpmag);
	  double z = (rmag - rpmag);

	  double U;
	  U = PairActionVector(PairIndex)->U(q,z,s2, level);
	  // Subtract off long-range part from short-range action
	  if (PairActionVector(PairIndex)->IsLongRange())
	    U -= 0.5* (PairActionVector(PairIndex)->Ulong(level)(rmag) +
		       PairActionVector(PairIndex)->Ulong(level)(rpmag));
	  TotalU += U;
	}
      }
    }
  }

  if (PathData.Path.LongRange){
    // Now add in the long-range part of the action
    // In primitive form, end slices get weighted by 1/2.  Others by 1.
    for (int slice=startSlice;slice<=endSlice;slice+=skip) 
      if ((slice==startSlice) || (slice == endSlice))
	TotalU += 0.5*LongRange_U (slice, level);
      else
	TotalU += LongRange_U (slice, level);
  }
  return (TotalU);

  if (PathData.Path.LongRange){
    // Now add in the long-range part of the action
    // In primitive form, end slices get weighted by 1/2.  Others by 1.
    if (UseRPA) {
      for (int slice=startSlice;slice<=endSlice;slice+=skip) 
	if ((slice==startSlice) || (slice == endSlice))
	  TotalU += 0.5*LongRange_U_RPA (slice, level);
	else
	  TotalU += LongRange_U_RPA (slice, level);
    }
    else {
      for (int slice=startSlice;slice<=endSlice;slice+=skip) 
	if ((slice==startSlice) || (slice == endSlice))
	  TotalU += 0.5*LongRange_U (slice, level);
	else
	  TotalU += LongRange_U (slice, level);
    }
  }
  return (TotalU);
}
