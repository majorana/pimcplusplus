#include "BisectionBlock.h"

void BisectionBlockClass::Read(IOSectionClass &in)
{
  string permuteType, speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  assert (in.ReadVar ("name", Name));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  IsFermion = PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION;

  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    permuteStage = new TablePermuteStageClass (PathData, SpeciesNum, NumLevels);
  else if (permuteType == "NONE") 
    permuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  permuteStage->Read (in);
  Stages.push_back (permuteStage);
  
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level);
    newStage->Actions.push_back(&PathData.Actions.ShortRange);
    //newStage->Actions.push_back(&PathData.Actions.LongRange);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    if ((level == 0) && (PathData.Actions.NodalActions(SpeciesNum)!=NULL))
      newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }

  // Add the second stage of the permutation step
  Stages.push_back (permuteStage);

  // Add the nodal action stage, if necessary
  
}


void BisectionBlockClass::ChooseTimeSlices()
{
  PathClass &Path = PathData.Path;
  int myProc = PathData.Communicator.MyProc();
  // do something special to avoid moving reference slice
  if (IsFermion &&
      Path.SliceOwner(Path.GetRefSlice()) == myProc) {
    
    int bSlices = 1<<NumLevels;
    int myStart, myEnd;
    Path.SliceRange (myProc, myStart, myEnd);
    // refSlice is relative to my first slice
    int refSlice = Path.GetRefSlice() - myStart;
    int numSlices = Path.NumTimeSlices();
    if (refSlice < bSlices) {
      int numStarts = numSlices - bSlices - refSlice;
      Slice1 = refSlice + Path.Random.LocalInt(numStarts);
      Slice2 = Slice1 + bSlices;
      assert (Slice2 < numSlices);
    }
    else if (refSlice > (numSlices-1-bSlices)) {
      int numStarts = refSlice - bSlices + 1;
      Slice1 = Path.Random.LocalInt(numStarts);
      Slice2 = Slice1+bSlices;
      assert (Slice2 <= refSlice);
    }
    else {
      int numBefore = refSlice - bSlices +1;
      int numAfter = numSlices -bSlices - refSlice;
      int numStarts = numBefore+numAfter;
      int start = Path.Random.LocalInt (numStarts);
      if (start < numBefore)
	Slice1 = start;
      else
	Slice1 = start-numBefore+refSlice;
      Slice2 = Slice1+bSlices;
      assert ((refSlice <= Slice1) || (refSlice >= Slice2));
      assert (Slice2 < numSlices);
    }
  }
  // Bryan should fill this part in.
  // else if (IsOpen(ptcl)) 
  else {
    int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    Slice1 = PathData.Path.Random.LocalInt (numLeft);
    Slice2 = Slice1+sliceSep;
  }
}

void BisectionBlockClass::MakeMove()
{
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0)=-1;
    MultiStageClass::MakeMove();
  }
}
