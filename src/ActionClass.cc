#include "ActionClass.h"

double ActionClass::calcTotalAction(Array<ParticleID,1> changedParticles,
				  int StartSlice, int EndSlice, int level)
{
  // First, sum the pair actions
  double TotalU = 0.0;
  double TotalK = 0.0;
  int NumChangedPtcls = changedParticles.size();
  int Species1, Species2, Ptcl1, Ptcl2;
  int NumSpecies;

  ArrayOfIdenticalParticlesClass &IdentPtcls = *myIdenticalParticleArray;
  NumSpecies = IdentPtcls.size();

  int skip = 1<<level;
  int levelTau = tau*1<<level;
  for (int i=0; i<NumChangedPtcls; i++)
    {
      Species1 = changedParticles(i)[0];
      Ptcl1 = changedParticles(i)[1];
      for (int Species2=0; Species2<NumSpecies; Species2++) {
	int NumPtcls2 = IdentPtcls(Species2).NumParticles;
	for (int Ptcl2=0; Ptcl2<NumPtcls2; Ptcl2++) {
	  for (int Slice=StartSlice; Slice < EndSlice; Slice+=skip) {	    
	    double NotMyself = (double)((Ptcl1!=Ptcl2)||(Species1!=Species2));
	    dVec r1 = IdentPtcls(Species1).Path(Ptcl1,Slice);
	    dVec r2 = IdentPtcls(Species2).Path(Ptcl2,Slice);
	    dVec rp1 = IdentPtcls(Species1).Path(Ptcl1,Slice+skip);
	    dVec rp2 = IdentPtcls(Species2).Path(Ptcl2,Slice+skip);
	    dVec r = r1 - r2;
	    dVec rp = rp1 - rp2;
	    double rmag = sqrt(dot(r,r));
	    double rpmag = sqrt(dot(rp,rp));
	    
	    double s = sqrt(dot (r-rp, r-rp));
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);
	    int PairIndex = PairMatrix(Species1, Species2);
	    TotalU += NotMyself*
	      PairActionVector(PairIndex).calcUrrptau(s,q,z, level);
	  }
	}
      }
      
      // Now, sum up the kinetic action
      double FourLambdaTauInv=1/(4*IdentPtcls(Species1).lambda*tau);
      for (int Slice=StartSlice; Slice < EndSlice; Slice+=skip) {
	dVec r1 = IdentPtcls(Species1).Path(Ptcl1,Slice);
	dVec r2 = IdentPtcls(Species1).Path(Ptcl1,Slice+skip);
	double LinkDistSqrd=distSqrd(r1,r2);  ///This function has to be written and possibly memoized or something?
	TotalK += LinkDistSqrd*FourLambdaTauInv; //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
	
	
      }
    
      return (TotalK+TotalU);
}
