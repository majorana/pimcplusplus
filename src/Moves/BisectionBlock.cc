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
    //newStage->Actions.push_back(&PathData.Actions.ShortRange);
    //newStage->Actions.push_back(&PathData.Actions.LongRange);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }
  // Add the second stage of the permutation step
  //  Stages.push_back (permuteStage);

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

    //  cerr<<"I have chosen "<<Slice1<<" and "<<Slice2<<endl;
  }

}

void BisectionBlockClass::MakeMove()
{
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);
  ActiveParticles.resize(1);
  ActiveParticles(0)=-1;
  for (int step=0; step<StepsPerBlock; step++)
    MultiStageLocalClass::MakeMove();
}
