#include "BisectionBlock.h"

void BisectionBlockClass::Read(IOSectionClass &in)
{
  string PermuteType;
  StageClass *permuteStage;
  assert (in.ReadVar ("PermuteType", PermuteType));
  if (PermuteType == "TABLE") {
    PermuteTableClass *newPermute = new PermuteTableClass (PathData);
    newPermute->Read (in);
    permuteStage = newPermute;
    Stages.push_back (newPermute);
  }
  else {
    cerr << "Unrecognized PermuteType, """ << PermuteType << """\n";
    exit(EXIT_FAILURE);
  }
  assert (in.ReadVar ("NumLevels", NumLevels));
  for (int level=NumLevels; level>=0; level--) {
    BisectionStageClass &newStage = new BisectionStageClass (pathData, level);
    Stages.push_back (newStage);
  }
  // Add the second stage of the permutation step
  Stages.push_back (permuteStage);

  // Add the nodal action stage, if necessary
  
}


