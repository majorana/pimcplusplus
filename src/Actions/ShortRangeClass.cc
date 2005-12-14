#include "ShortRangeClass.h"
#include "../PathDataClass.h"
#include "ShortRangeOnClass.h"
#include "time.h"
#include <Common/MatrixOps/MatrixOps.h>

///This has to be called after pathdata knows how many
///particles it has
void ShortRangeClass::Read(IOSectionClass& in)
{
  DoPtcl.resize(PathData.Path.NumParticles());
  TotalTime=0;
}

ShortRangeClass::ShortRangeClass(PathDataClass &pathData,
				 Array<PairActionFitClass* ,2> &pairMatrix) : 
  ToCheck(pathData,pairMatrix),
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix),
  m(2), NumBasisFuncs(4), Router(0.6), UseLowVariance(/*true*/false)
{
  Setup_ck();
}

void ShortRangeClass::Setup_ck()
{
  ck.resize(NumBasisFuncs);
  Array<double,1> h(NumBasisFuncs);
  Array<double,2> S(NumBasisFuncs);
  double Rtojplus1 = Router*Router;
  for (int j=1; j<=NumBasisFuncs; j++) {
    h(j-1) = Rtojplus1/(1.0+j);
    Rtojplus1 *= Router;
    for (int k=1; k<=NumBasisFuncs; k++)
      S(k-1,j-1) = pow(Router, m+k+j+1)/(m+k+j+1.0);
  }
  GJInverse (S);
  ck = 0.0;
  for (int i=0; i < NumBasisFuncs; i++)
    for (int j=0; j < NumBasisFuncs; j++)
      ck(i) += S(i,j) * h(j);

  /// HACK HACK HACK
  FILE *fout = fopen ("gr.dat", "w");
  for (double r=0.0; r<=Router; r+=0.001)
    fprintf (fout, "%1.18e %1.18e\n", r, g(r));
  fclose (fout);

}

inline double
ShortRangeClass::g(double r)
{
  if (r > Router)
    return 1.0;
  double rtokplusm = r;
  for (int i=0; i<m; i++)
    rtokplusm *= r;

  double gval = 0.0;
  for (int k=1; k<=NumBasisFuncs; k++) {
    gval += ck(k-1) * rtokplusm;
    rtokplusm *= r;
  }
  return gval;
}

double 
ShortRangeClass::SingleAction (int slice1, int slice2,
			       const Array<int,1> &changedParticles,
			       int level)
{
  double TotalU=0.0;
  //  int startTime=clock();
  //  for (int toRun=0;toRun<1000;toRun++){
  PathClass &Path=PathData.Path;
  // First, sum the pair actions
  for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
    Path.DoPtcl(ptcl)=true;
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
// 	  if (isinf(U)) {
// 	    cerr << "inf at slice=" << slice << " ptcl1=" << ptcl1
// 		 << " ptcl2=" << ptcl2 << endl;
// 	    cerr << "r = " << r << endl;
// 	    cerr << "rp = " << rp << endl;
// 	    cerr << "r1 = " << PathData.Path(slice+skip, ptcl1) << endl;
// 	    cerr << "r2 = " << PathData.Path(slice+skip, ptcl2) << endl;
// 	  }

	  //	  if (ptcl2==4 && ptcl1==0 && slice==8)
	    //	    cerr<<"MY real U is "<<U<<endl;
	  //	  if (isnan(U)){
	  //	    cerr << "Before long range sub:  ptcl1=" << ptcl1
	  //		 << " ptcl2=" << ptcl2 << " slice="<< slice 
	  //		 << "q: "<<q<<"z: "<<z<<"s2: "<<s2<<"level: "<<level<<endl;
	  //	    cerr<<"U is "<<U<<endl;
	  //	  }
	  // Subtract off long-range part from short-range action
	  if (PA.IsLongRange() && PathData.Actions.UseLongRange)
	    U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));

// 	  if ((slice1==0) && (slice2 == Path.NumTimeSlices()))
// 	      cerr << "U = " << U << endl;
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
//  cerr<<endl;
  //  cerr<<"My U is "<<TotalU<<endl;
  //  cerr<<"Other U is "<<ToCheck.Action(slice1, slice2,
  //				      changedParticles,
  //				      level);
  //cerr<<endl;
  //cerr<<endl;
  //  }
  //  int endTime=clock();
//   cerr<<"Diff is "<< (double)(endTime-startTime)/(double)CLOCKS_PER_SEC<<endl;
//   startTime=clock();
//   double dummyVar;
//   for (int numTime=0;numTime<1000;numTime++)
//     // assert(TotalU-ToCheck.Action(slice1, slice2,changedParticles,level)<1e-2);
//     dummyVar+=ToCheck.Action(slice1, slice2,changedParticles,level);
//   endTime=clock();
//   cerr<<"O(n) Diff is "<<(double)(endTime-startTime)/(double)CLOCKS_PER_SEC<<endl;
//   cerr<<endl;
//   cerr<<"A: "<<dummyVar<<" "<<TotalU<<endl;
  return (TotalU);
}



double 
ShortRangeClass::d_dBeta (int slice1, int slice2, int level)
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
	if (pa.IsLongRange() && PathData.Actions.UseLongRange)
	  dU -= 0.5*(pa.dUlong(level)(rmag)+pa.dUlong(level)(rpmag));
      }
    }
  }
  return dU;
}


void
ShortRangeClass::GradAction(int slice1, int slice2, 
			    const Array<int,1> &ptcls, int level,
			    Array<dVec,1> &gradVec)
{
  PathClass &Path = PathData.Path;
  int skip = (1<<level);
  assert (gradVec.size() == ptcls.size());
  for (int pi=0; pi<ptcls.size(); pi++) {
    int ptcl1 = ptcls(pi);
    int species1 = Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      PairActionFitClass &PA=*(PathData.Actions.PairMatrix(species1,species2));
      if (ptcl1 != ptcl2) {
	for (int slice=slice1; slice<slice2; slice += skip) {
	  	  dVec r, rp;
	  double rmag, rpmag, du_dq, du_dz;
	  Path.DistDisp(slice, slice+skip, ptcl1, ptcl2, rmag, rpmag, r, rp);
	  double q = 0.5*(rmag+rpmag);
	  double z = (rmag-rpmag);
	  double s2 = dot (r-rp, r-rp);
	  PA.Derivs(q,z,s2,level,du_dq, du_dz);
	  Vec3 rhat  = (1.0/rmag)*r;
	  Vec3 rphat = (1.0/rpmag)*rp;
	  
	  double g1 = 1.0;
	  double g2 = 1.0;
	  if (UseLowVariance) {
	    g1 = g(rmag);
	    g2 = g(rpmag);
	  }
	  gradVec(pi) -= (g1*(0.5*du_dq + du_dz)*rhat + 
			  g2*(0.5*du_dq - du_dz)*rphat);
	  // gradVec(pi) -= (0.5*du_dq*(rhat+rphat) + du_dz*(rhat-rphat));
	  /// Now, subtract off long-range part that shouldn't be in
	  /// here 
	  if (PA.IsLongRange() && PathData.Actions.UseLongRange)
	    gradVec(pi) += 0.5*(PA.Ulong(level).Deriv(rmag)*g1*rhat+
				PA.Ulong(level).Deriv(rpmag)*g2*rphat);
	}
      }
    }
  }
}
