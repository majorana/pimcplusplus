#include "BlockMove.h"

double CycleBlockMoveClass::AcceptanceRatio()
{
  return (double)NumAccepted/(double)NumMoves;
}

void CycleBlockMoveClass::MakeMove()
{
  // First, decide on the chunk of slices we're working on
  int slice1;
  int slice2;
  int numSlices=PathData.NumTimeSlices();

  if (!PathData.Path.OpenPaths){
    double xi=PathData.Path.Random.Local();
    slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
    slice1=slice2-(1<<NumLevels);
  }
  else {
    // First, decide on the chunk of slices we're working on
    do{
      double xi=PathData.Path.Random.Local();
      slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
      slice1=slice2-(1<<NumLevels);
      //      if (slice2>=PathData.Path.NumTimeSlices()-1){
	//	cerr<<"ERROR!"<<slice2<<endl;
      //	slice2--;
      //      }
    } while ((slice1<=(int)(PathData.Path.OpenLink) && (int)(PathData.Path.OpenLink)<=slice2));
  }
  //  if (slice2>=PathData.Path.NumTimeSlices()-1){
  //    cerr<<"ERROR! ERROR! ERROR!";
  //  }
  //  cerr<<slice1<<" "<<slice2<<" "<<PathData.Path.NumTimeSlices()<<endl;
  PathData.MoveJoin(slice2);
  
  int step = 0;
  // Now, construct the Forward table
  Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);
  
  int NumPerms = 0;
  for (int step=0; step<StepsPerBlock; step++) {
    // Choose a permutation cycle
    double forwT = Forw->AttemptPermutation();
    double revT = Rev->CalcReverseProb(*Forw);
    double Tratio = forwT/revT;
    int len=Forw->CurrentCycle.Length;
    // double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1];
    double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
    double psi = PathData.Path.Random.Local();
    Array<int,1> currentParticles=Forw->CurrentParticles();  
    double pi_ratio = exp(-actionChange);
    double accept_prob = min(1.0, pi_ratio/Tratio);
    //    cerr<<"My accept prob is "<<pi_ratio/Tratio<<endl;
    //    if (log(psi) < -(actionChange + log(Tratio))) {
    if  (accept_prob > psi) {
      bool acceptBisect = 
	Bisection.Bisect(slice1, NumLevels, currentParticles, actionChange);
      if (acceptBisect){
	PathData.AcceptMove(slice1,slice2,currentParticles);

	// We don't construct a new table for single-ptcl moves!
	if (Forw->CurrentCycle.Length!=1){
	  NumPerms++;
	  PermuteTableClass* tempPtr=Forw;
	  Forw=Rev;
	  Rev=tempPtr;
	  cerr<<"I've accepted a move of "<<Forw->CurrentCycle.Length<<" size"<<endl;

	}

	NumAccepted++;
      }
      else{
	PathData.RejectMove(slice1,slice2,currentParticles);
      }
    }
    else
      PathData.RejectMove(slice1, slice2,currentParticles);
  }
  NumMoves+=StepsPerBlock;
}

void CycleBlockMoveClass::Read(IOSectionClass  &in)
{
  Forw->Read(in);
  Rev->Read(in);
  string speciesName;
  assert(in.ReadVar("Species",speciesName));
  assert(in.ReadVar("NumLevels",NumLevels));
  assert(in.ReadVar("Steps",StepsPerBlock));
  SpeciesNum=PathData.SpeciesNum(speciesName);
  assert(in.ReadVar("name",Name));
}

