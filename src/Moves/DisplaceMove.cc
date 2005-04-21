#include "DisplaceMove.h"

double DisplaceStageClass::Sample (int &slice1, int &slice2,
				   Array<int,1> &activeParticles)
{
  /// Now, choose a random displacement 
  dVec disp;
  PathData.Path.Random.CommonGaussianVec (Sigma, disp);

  // Actually displace the path
  SetMode(NEWMODE);
  int ptcl = activeParticles(0);
  for (int slice=0; slice<PathData.NumTimeSlices(); slice++)
    PathData.Path(slice, ptcl) = PathData.Path(slice, ptcl) + disp;

  // And return sample probability ratio
  return 1.0;
}

void
DisplaceMoveClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar ("Sigma", Sigma));
  DisplaceStage.Sigma = Sigma;
  assert(in.ReadVar("name",Name));
  Array<string,1> activeSpeciesNames;

  // Read in the active species.
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  Array<int,1> activeSpecies(activeSpeciesNames.size());
  for (int i=0; i<activeSpecies.size(); i++)
    activeSpecies(i) = PathData.Path.SpeciesNum(activeSpeciesNames(i));
  SetActiveSpecies (activeSpecies);

  // Only move 1 particle at a time.
  SetNumParticlesToMove(1);

  // Construct action list
  DisplaceStage.Actions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange) 
    if (PathData.Actions.UseRPA)
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    else
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRange);

//   for (int i=0; i<PathData.Actions.NodalActions.size(); i++)
//     DisplaceStage.Actions.push_back(PathData.Actions.NodalActions(i));
  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      cerr << "Adding fermion node action for species " 
	   << activeSpeciesNames(i) << endl;
      DisplaceStage.Actions.push_back
	(PathData.Actions.NodalActions(speciesNum));
    }
  }
    
  // Now construct stage list
  Stages.push_back(&DisplaceStage);
}

void
DisplaceMoveClass::MakeMove ()
{
  // Move the Join out of the way.
  PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);

  // First, choose particle to move
  int numActive = 0;
  for (int i=0; i<ActiveSpecies.size(); i++)
    numActive += PathData.Path.Species(ActiveSpecies(i)).NumParticles;
  ActiveParticles.resize(1);


  int index = PathData.Path.Random.CommonInt(numActive);
  bool done = false;
  int speciesIndex = 0;
  while ((speciesIndex < ActiveSpecies.size()) && (!done)) {
    SpeciesClass &species = PathData.Path.Species(ActiveSpecies(speciesIndex));
    if (index < species.NumParticles) {
      ActiveParticles(0) = index+species.FirstPtcl;
      done = true;
    }
    else
      index -= species.NumParticles;
    speciesIndex++;
  }

  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;

  // Now call MultiStageClass' MakeMove
  MultiStageClass::MakeMove();
}

void
DisplaceMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();

}
