#include "BisectionBlock.h"

void BisectionBlockClass::Read(IOSectionClass &in)
{
  string permuteType, speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);


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
  
  for (int level=NumLevels; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level);
    Stages.push_back (newStage);
  }
  // Add the second stage of the permutation step
  Stages.push_back (permuteStage);

  // Add the nodal action stage, if necessary
  
}


void BisectionBlockClass::ChooseTimeSlices()
{
  if (IsFermion) {
    // do something special to avoid moving reference slice
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
  for (int step=0; step<StepsPerBlock; step++)
    MultiStageLocalClass::MakeMove();
}
