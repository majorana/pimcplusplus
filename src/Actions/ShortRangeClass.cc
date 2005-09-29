#include "ShortRangeClass.h"
#include "../PathDataClass.h"

///This has to be called after pathdata knows how many
///particles it has
void ShortRangeClass::Read(IOSectionClass& in)
{
  DoPtcl.resize(PathData.Path.NumParticles());
}

ShortRangeClass::ShortRangeClass(PathDataClass &pathData,
				 Array<PairActionFitClass* ,2> &pairMatrix) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix)
{
}

double ShortRangeClass::Action (int slice1, int slice2,
				const Array<int,1> &changedParticles,
				int level)
{
  PathClass &Path=PathData.Path;
  // First, sum the pair actions
  for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
    Path.DoPtcl(ptcl)=true;
  double TotalU = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = changedParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);

    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++) {
      if (Path.DoPtcl(ptcl2)){
	int species2=Path.ParticleSpeciesNum(ptcl2);
	PairActionFitClass &PA = *(PairMatrix(species1, species2));
	for (int slice=slice1;slice<slice2;slice+=skip){
	  dVec r, rp;
	  double rmag, rpmag;

	  PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				 rmag, rpmag, r, rp);
	  double s2 = dot (r-rp, r-rp);
	  double q = 0.5 * (rmag + rpmag);
	  double z = (rmag - rpmag);
	  double U;
	  U = PA.U(q,z,s2, level);
	  //	  if (ptcl2==4 && ptcl1==0 && slice==8)
	    //	    cerr<<"MY real U is "<<U<<endl;
	  //	  if (isnan(U)){
	  //	    cerr << "Before long range sub:  ptcl1=" << ptcl1
	  //		 << " ptcl2=" << ptcl2 << " slice="<< slice 
	  //		 << "q: "<<q<<"z: "<<z<<"s2: "<<s2<<"level: "<<level<<endl;
	  //	    cerr<<"U is "<<U<<endl;
	  //	  }
	  // Subtract off long-range part from short-range action
	  if (PA.IsLongRange())
	    U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
	  //	  if (isnan(U))
	    //	    cerr << "After  long range sub:  ptcl1=" << ptcl1
	    //		 << " ptcl2=" << ptcl2 << " slice="<< slice << endl;

	   // Code for vancacy project commented out
// 	  if (PathData.Path.ExistsCoupling>=-0.05 &&
// 	      (species1==1 || species2==1))
// 	    TotalU += sqrt(PathData.Path.ExistsCoupling)*U;
// 	  else
// 	    TotalU += U;

	  TotalU += U;
	}
      }
    }
  }

//   if (PathData.Path.LongRange){
//     // Now add in the long-range part of the action
//     // In primitive form, end slices get weighted by 1/2.  Others by 1.
//     for (int slice=slice1;slice<=slice2;slice+=skip) 
//       if ((slice==slice1) || (slice == slice2))
// 	TotalU += 0.5*LongRange_U (slice, level);
//       else
// 	TotalU += LongRange_U (slice, level);
//   }
//   return (TotalU);

//   if (PathData.Path.LongRange){
//     // Now add in the long-range part of the action
//     // In primitive form, end slices get weighted by 1/2.  Others by 1.
//     if (UseRPA) {
//       for (int slice=slice1;slice<=slice2;slice+=skip) 
// 	if ((slice==slice1) || (slice == slice2))
// 	  TotalU += 0.5*LongRange_U_RPA (slice, level);
// 	else
// 	  TotalU += LongRange_U_RPA (slice, level);
//     }
//     else {
//       for (int slice=slice1;slice<=slice2;slice+=skip) 
// 	if ((slice==slice1) || (slice == slice2))
// 	  TotalU += 0.5*LongRange_U (slice, level);
// 	else
// 	  TotalU += LongRange_U (slice, level);
//     }
//   }
//  for (int counter=0;counter<TotalUArray.size();counter++){
//    cerr<<"My link "<<counter+1<<" "<<counter+2<<" is "
//	<<TotalUArray(counter)<<endl;
//  }
  return (TotalU);
}



double ShortRangeClass::d_dBeta (int slice1, int slice2,
				 int level)
{
  double levelTau=Path.tau;
  int skip = 1<<level;
  //  int slice2 = slice1 + (1<<level);
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  double dU = 0.0;
  for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      for (int slice=slice1;slice<slice2;slice+=skip){
	dVec r, rp;
	double rmag, rpmag;
	PathData.Path.DistDisp(slice,slice+skip,ptcl1,ptcl2,rmag,rpmag,r,rp);
	
	double s2 = dot(r-rp, r-rp);
	double q = 0.5*(rmag+rpmag);
	double z = (rmag-rpmag);
	
	PairActionFitClass& pa=
	  *(PairMatrix(species1, PathData.Path.ParticleSpeciesNum(ptcl2)));
	dU += pa.dU(q, z, s2, level);
	// Subtract off long-range part from short-range action
	if (pa.IsLongRange())
	  dU -= 0.5*(pa.dUlong(level)(rmag)+pa.dUlong(level)(rpmag));
      }
    }
  }
  return dU;
}
