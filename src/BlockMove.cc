#include "BlockMove.h"
double CycleBlockMoveClass::AcceptanceRatio()
{
  return (double)NumAccepted/(double)NumMoves;
}

void CycleBlockMoveClass::MakeMove()
{
  // First, decide on the chunk of slices we're working on
  int numSlices=PathData.NumTimeSlices();
  double xi=PathData.Path.Random.Local();
  int slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
  int slice1=slice2-(1<<NumLevels);
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
    double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
    double psi = PathData.Path.Random.Local();
    Array<int,1> currentParticles=Forw->CurrentParticles();  
    if (log(psi) < (-actionChange - log(Tratio))) {
      bool acceptBisect = 
	Bisection.Bisect(slice1, NumLevels, currentParticles, actionChange);
      if (acceptBisect){
	PathData.AcceptMove(slice1,slice2,currentParticles);

	// We don't construct a new table for single-ptcl moves!
	//Right now we are!
// 	if (Forw->CurrentCycle.Length!=1){
// 	  NumPerms++;
// 	  PermuteTableClass* tempPtr=Forw;
// 	  Forw=Rev;
// 	  Rev=tempPtr;
// 	}
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

