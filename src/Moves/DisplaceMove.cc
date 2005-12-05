#include "DisplaceMove.h"

double DisplaceStageClass::Sample (int &slice1, int &slice2,
				   Array<int,1> &activeParticles)
{
  /// Now, choose a random displacement 
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    dVec disp;
    PathData.Path.Random.CommonGaussianVec (Sigma, disp);

    // Actually displace the path
    SetMode(NEWMODE);
    for (int slice=0; slice<PathData.NumTimeSlices(); slice++)
      PathData.Path(slice, ptcl) = PathData.Path(slice, ptcl) + disp;
  }

  // And return sample probability ratio
  return 1.0;
}

void
DisplaceMoveClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar ("Sigma", Sigma));
  DisplaceStage.Sigma = Sigma;
  Array<string,1> activeSpeciesNames;

  int numToMove;
  assert(in.ReadVar("NumToMove", numToMove));
  SetNumParticlesToMove(numToMove);

  // Read in the active species.
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  Array<int,1> activeSpecies(activeSpeciesNames.size());
  for (int i=0; i<activeSpecies.size(); i++)
    activeSpecies(i) = PathData.Path.SpeciesNum(activeSpeciesNames(i));
  SetActiveSpecies (activeSpecies);

//   // Move all particles at the same time.
//   int totalNum = 0;
//   for (int si=0; si<activeSpecies.size(); si++)
//     totalNum += PathData.Path.Species(activeSpecies(si)).NumParticles;
//   SetNumParticlesToMove(totalNum);

  // Construct action list
  DisplaceStage.Actions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange) 
    if (PathData.Actions.UseRPA)
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    else
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRange);

  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      cerr << "DisplaceMove adding fermion node action for species " 
	   << activeSpeciesNames(i) << endl;
      DisplaceStage.Actions.push_back
	(PathData.Actions.NodalActions(speciesNum));
    }
  }
    
  // Now construct stage list
  Stages.push_back(&DisplaceStage);

  ActiveParticles.resize(NumParticlesToMove);
}

void
DisplaceMoveClass::MakeMove ()
{
  // Move the Join out of the way.
  PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);

  vector<int> ptclList;
  int numLeft=0;
  for (int i=0; i<ActiveSpecies.size(); i++) {
    SpeciesClass &species = PathData.Path.Species(ActiveSpecies(i));
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      ptclList.push_back(ptcl);
      numLeft++;
    }
  }
  
  // First, choose particle to move
  for (int i=0; i<NumParticlesToMove; i++) {
    int index = PathData.Path.Random.CommonInt(numLeft);
    vector<int>::iterator iter = ptclList.begin();
    ActiveParticles(i) = ptclList[index];
    for (int j=0; j<index; j++)
      iter++;
    ptclList.erase(iter);
    numLeft--;
  }
  cerr << "ActiveParticles = " << ActiveParticles << endl;

//   int ptclIndex = 0;
//   for (int si=0; si < ActiveSpecies.size(); si++) {
//     SpeciesClass &species = PathData.Path.Species(ActiveSpecies(si));
//     for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
//       ActiveParticles(ptclIndex) = ptcl;
//       ptclIndex++;
//     }
//   }
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;

  // Now call MultiStageClass' MakeMove
  MultiStageClass::MakeMove();
}
