#include "PermuteTableClass.h"


// void CycleClass::CanonicalPermRep(Array<int,1> Perm)
// {
//   // Set to identity permutation
//   for (int i=0;i<myArray;i++)
//     Perm(i)=i;
//   // Apply cyclic permutation
//   for(int i=0;i<Ncycles;i++)
//     Perm(CycleRep(i)) = CycleRep((i+1) % Ncycles);
// }
    


void PermuteTableClass::ConstructHTable()
{
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  double lambda = PathData.Species(SpeciesNum).lambda;
  double beta = PathData.Action.tau * (double) (Slice2-Slice1);
  double fourLambdaBetaInv = 1.0/(4.0*lambda*beta);
  int N = lastPtcl-firstPtcl+1;
  HTable.resize(N,N);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      dVec disp_ij = PathData(Slice2,j+firstPtcl)-PathData(Slice1,i+firstPtcl);
      PathData.DistanceTable->PutInBox(disp_ij);
      double dist_ij = dot (disp_ij,disp_ij);
      dVec disp_ii = PathData(Slice2,i+firstPtcl)-PathData(Slice1,i+firstPtcl);
      PathData.DistanceTable->PutInBox(disp_ii);
      double dist_ii = dot (disp_ii,disp_ii);
      HTable(i,j) = exp((-dist_ij + dist_ii)*fourLambdaBetaInv);
      if (HTable(i,j) < epsilon)
	HTable(i,j) = 0.0;
    }
    // Now we should really sort with respect to j, but we won't for
    // now because we are far too lazy and it's a Friday.
  }
}

 
///// This constructs the Htable for the reverse move from the Htable
///// for the forward move.  i.e. Constructs Hrev(i,j) = H(P(i),P(j)).
// void PermuteTableClass::PermuteHTable(const CycleClass &myPerm,
// 				      const PermuteTableClass &forwardTable)
// {
//   int N=HTable.extent(0);
//   Array<int,1> P(PathData.NumParticles());
//   myPerm.CanonicalPermRep(P);
//   for (int i=0; i<N; i++) {
//     int P_i = P(i);
//     for (int j=0; j<N; j++) {
//       int P_j = P(j);
//       HTable(P_i,P_j) = forwardTable.HTable(i,j);
//     }
//   }
 
// } 
 
 
// Actually applies my cyclic permutation to a path at a given time
// slice and the permutation vector in the path
void CycleClass::Apply(PathClass &path, int firstPtcl, int slice)
{
  dVec tempPos = path(slice, CycleRep(0)+firstPtcl);
  int tempPtcl = path.Permutation(CycleRep(0)+firstPtcl);
  for(int i=0;i<Ncycles-1;i++) {
    path.SetPos(slice, CycleRep(i)+firstPtcl, 
		path(slice,CycleRep(i+1)+firstPtcl));
    path.Permutation.Set(CycleRep(i)+firstPtcl,
			 path.Permutation(CycleRep(i+1)+firstPtcl));
  }
  path.SetPos(slice,CycleRep(Ncycles-1)+firstPtcl,tempPos);
  path.Permutation.Set(CycleRep(Ncycles-1)+firstPtcl,tempPtcl);
}


double PermuteTableClass::AttemptPermutation()
{
  //Get a random number number from the local processor stream.
  
  double xi=PathData.Path.Random.Local(); 
  int index=FindEntry(xi);
  CycleClass &cycle = PermTable(index);
  
  // Now, apply the permutation to the Path
  int firstPtcl=PathData.Species(SpeciesNum).FirstPtcl;
  cycle.Apply(PathData.Path,firstPtcl,Slice2);
  return (cycle.P * NormInv);
}

double PermuteTableClass::CalcReverseProb(const CycleClass &myPerm,
				     const PermuteTableClass &forwardTable)
{

  //We reconstruct things from scratch to make sure we do the right
  //thing. We can try to incrementally make this faster.
  ConstructCycleTable(forwardTable.SpeciesNum,
		      forwardTable.Slice1,forwardTable.Slice2); 
  return 1/(myPerm.P)*NormInv;
	 
	   
}


void PermuteTableClass::Read(IOSectionClass &inSection)
{
  Array<double,1> tempGamma;
  assert(inSection.ReadVar("Gamma",tempGamma));
  assert(tempGamma.size() == 4);
  for (int counter=0;counter<tempGamma.size();counter++){
    Gamma(counter)=tempGamma(counter);
  }
  assert(inSection.ReadVar("epsilon",epsilon));
  //  assert(inSection.ReadVar("SpeciesNum",SpeciesNum));
  
  
}


void PermuteTableClass::ConstructCycleTable(int speciesNum,
					    int slice1, int slice2)
{
  Slice1=slice1;
  Slice2=slice2;
  SpeciesNum=speciesNum;
  ConstructHTable();
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  PermTable.resize(TableSize);
  NumEntries = 0;
  //  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  //  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  int N=HTable.extent(0);
  for (int i=0; i<N; i++) {
    ///Single particle move
    CycleClass tempPerm;
    double hprod=1.0;
    tempPerm.Ncycles=1;
    tempPerm.CycleRep[0]=i;
    tempPerm.P=Gamma(0);
    tempPerm.C=tempPerm.P;
    if (NumEntries!=0){
      tempPerm.C+=PermTable(NumEntries-1).C;
    }
    AddEntry(tempPerm);
    for (int j=i+1; j<N; j++) {// 2 and higher cycles
      //2 cycle permutations
      tempPerm.CycleRep[1]=j;
      hprod*=HTable(i,j);
      tempPerm.Ncycles=2;
      tempPerm.P=Gamma(1)*hprod*HTable(j,i);
      tempPerm.C=tempPerm.P+PermTable(NumEntries-1).C;      
      if (hprod != 0.0) {
	AddEntry(tempPerm);
	for (int k=i+1;k<N;k++){//3 and higher cycles
	  //3 cycle permutations
	  if (k!=j){
	    tempPerm.CycleRep[2]=k;
	    hprod*=HTable(j,k);
	    tempPerm.Ncycles=3;
	    tempPerm.P=Gamma(2)*hprod*HTable(k,i);
	    tempPerm.C=tempPerm.P+PermTable(NumEntries-1).C;
	    if (hprod != 0.0) {
	      AddEntry(tempPerm);
	      for (int l=i+1;l<N;l++)
		if ((l!=j) && (l!=k)){
		  hprod*=HTable(k,l);
		  tempPerm.CycleRep[3]=l;
		  tempPerm.Ncycles=4;
		  tempPerm.P=Gamma(3)*hprod*HTable(l,i);
		  tempPerm.C=tempPerm.P+PermTable(NumEntries-1).C;
		  if (tempPerm.P != 0.0) 
		    AddEntry(tempPerm);
		}
	    }
	  }
	}
      }
    }
  }
  Norm = PermTable(NumEntries-1).C;
  NormInv = 1.0/Norm;
}
