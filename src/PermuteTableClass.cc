#include "PermuteTableClass.h"

void PermuteTableClass::ConstructHTable(int slice)
{
  int Slice2 = Slice1 + (1<<NumLevels);
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  double lambda = PathData.Species(SpeciesNum).lambda;
  double beta = PathData.Action.tau * (double) (1<<NumLevels);
  double fourLambaBetaInv = 1.0/(4.0*lambda*beta);
  HTable.resize(lastPtcl-firstPtcl+1,lastPtcl-firstPtcl+1);
  for (int i=firstPtcl; i<=lastPtcl; i++) {
    for (int j=firstPtcl; <=lastPtcl; j++) {
      dVec disp_ij = PathData(Slice2,j)-PathData(Slice1,i);
      PathData.DistanceTable->PutInBox(disp_ij);
      double dist_ij = dot (disp_ij,disp_ij);
      dVec disp_ii = PathData(Slice2,i)-PathData(Slice1,i);
      PathData.DistanceTable->PutInBox(disp_ii);
      double dist_ii = dot (disp_ii,disp_ii);
      HTable(i,j).j=j;
      HTable(i,j)=exp((-dist_ij + dist_ii)*fourLambdaBetaInv);
      if (HTable(i,j) < epsilon)
	HTable(i,j) = 0.0;
    }
    // Now we should really sort with respect to j, but we won't for
    // now because we are far too lazy and it's a Friday.
  }
}




void PermuteTableClass::CalcPermProbs()
{

  PermProbs.resize(tableSize);
  NumEntries = 0;
  int firstPtcl = PathData.Species(SpeciesNum).FirstPtcl;
  int lastPtcl = PathData.Species(SpeciesNum).LastPtcl;
  for (int i=firstPtcl; i<=lastPtcl; i++) {
    ///Single particle move
    PermClass tempPerm;
    double hprod=1.0;
    tempPerm.Ncycles=1;
    tempPerm.Perm[0]=i;
    tempPerm.P=Gamma[0];
    tempPerm.C=tempPerm.P;
    if (NumEntries!=0){
      tempPerm.P+=PermTable(NumEntries-1).C;
    }
    AddEntry(temPerm);
    for (int j=i+1; j<=lastPtcl; j++) {// 2 and higher cycles
      //2 cycle permutations
      tempPerm.Perm[1]=j;
      hprod*=HTable(i,j);
      tempPerm.Ncycles=2;
      tempPerm.P=gamma[1]*hprod*HTable(j,i);
      tempPerm.C=tempPerm.P+PermTable(NumEntries-1).C;      
      if (hprod != 0.0) {
	AddEntry(tempPerm);
	for (int k=i+1;k<=lastPtcl;k++){//3 and higher cycles
	  //3 cycle permutations
	  if (k!=j){
	    tempPerm.Perm[2]=k;
	    hprod*=HTable(j,k);
	    tempPerm.Ncycles=3;
	    tempPerm.P=gamma[2]*hprod*HTable(k,i);
	    tempPerm.C=tempPerm.P+PermTable(NumEntries-1).C;
	    if (hprod != 0.0) {
	      AddEntry(tempPerm);
	      for (l=i+1;l<=lastPtcl;l++){
		if ((l!=j) && (l!=k)){
		  hprod*=HTable(k,l);
		  tempPerm.Perm[3]=l;
		  tempPerm.Ncycles=4;
		  tempPerm.P=gamma[3]*hProd*HTable(l,i);
		  if (tempPerm.P != 0.0) 
		    AddEntry(tempPerm);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  double norm = 1.0/PermTable(NumEntries-1).C;
  for (int i=0; i<NumEntries; i++) {
    PermTable(i).P *= norm;
    PermTable(i).C *= norm;
  }
}
