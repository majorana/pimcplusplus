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
  int slice2=(int)(xi*(double)(numSlices-(1<<NumLevels))+(1<<NumLevels));
  int slice1=slice2-(1<<NumLevels);
  PathData.MoveJoin(slice2);
  
  int step = 0;
  // Now, construct the Forward table
  Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);

  for (int step=0; step<StepsPerBlock; step++) {
    //    cerr<<"I'm looping"<<endl;
    // Choose a permutation cycle
    double forwProb = Forw->AttemptPermutation();
    Array<int,1> currentParticles=Forw->CurrentParticles();
    bool acceptBisect = Bisection.Bisect(slice1,NumLevels,
					 currentParticles);
    if (acceptBisect){
      //      cerr<<"I'm accepting!"<<endl;
      double revProb=Rev->CalcReverseProb(*Forw);
      if (revProb/forwProb > PathData.Path.Random.Local()){
	PathData.AcceptMove(slice1,slice2,currentParticles);
	// We don't construct a new table for single-ptcl moves!
	if (Forw->CurrentCycle.Length!=1){
	  PermuteTableClass* tempPtr=Forw;
	  Forw=Rev;
	  Rev=tempPtr;
	}
	NumAccepted++;
      }
      else {
	PathData.RejectMove(slice1,slice2,currentParticles);
      }
    }
    else{
      PathData.RejectMove(slice1,slice2,currentParticles);
    }
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

