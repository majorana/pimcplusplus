#include "IonDisplaceMove.h"


double IonDisplaceStageClass::Sample (int &slice1, int &slice2,
				      Array<int,1> &activeParticles)
{
  SpeciesClass &ionSpecies = Path.Species(IonSpeciesNum);
  int first = ionSpecies.FirstPtcl;
  int last  = ionSpecies.LastPtcl;
  int N     = last - first + 1;
  if (DeltaRions.size() != N) {
    DeltaRions.resize(N);
    Weights.resize(N);
  }
  dVec zero(0.0);
  DeltaRions = zero;
  /// Now, choose a random displacement 
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    PathData.Path.Random.CommonGaussianVec (Sigma, DeltaRions(ptcl-first));

    // Actually displace the path
    SetMode(NEWMODE);
    for (int slice=0; slice<PathData.NumTimeSlices(); slice++)
      PathData.Path(slice, ptcl) = 
	PathData.Path(slice, ptcl) + DeltaRions(ptcl-first);
  }

  if (WarpElectrons)
    return DoElectronWarp();
  else
    return 1.0;

  // And return sample probability ratio
  return 1.0;
}

double
IonDisplaceStageClass::DoElectronWarp()
{
  SpeciesClass &ionSpecies = Path.Species(IonSpeciesNum);
  int ionFirst = ionSpecies.FirstPtcl;
  int ionLast  = ionSpecies.LastPtcl;
  int N     = ionLast - ionFirst + 1;
  dVec disp;
  double dist;
  for (int si=0; si<Path.NumSpecies(); si++) {
    SpeciesClass &species = Path.Species(si);
    if (species.lambda != 0.0) {
      for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
	for (int elec=species.FirstPtcl; elec<=species.LastPtcl; elec++) {
	  SetMode (OLDMODE);
	  Weights = 0.0;
	  double totalWeight = 0.0;
	  for (int ion=ionFirst; ion <= ionLast; ion++) {
	    Path.DistDisp(slice, elec, ion, dist, disp);
	    Weights(ion-ionFirst) = 1.0/(dist*dist*dist*dist);
	    totalWeight *= Weights(ion-ionFirst);
	  }
	  Weights *= (1.0/totalWeight);
	  SetMode (NEWMODE);
	  for (int ion=ionFirst; ion <= ionLast; ion++) 
	    Path(slice, elec) += Weights(ion-ionFirst)*DeltaRions(ion-ionFirst);
	}
      }
    }
  }
  return 1.0;
}


void
IonDisplaceMoveClass::MakeMove ()
{
  // Move the Join out of the way.
  PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);

  // First, displace the particles in the new copy
  SetMode(NEWMODE);

  // Construct list of particles to move
  vector<int> ptclList;
  int numLeft=0;
  SpeciesClass &species = PathData.Path.Species(IonSpeciesNum);
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    ptclList.push_back(ptcl);
    numLeft++;
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
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;
  // Now call MultiStageClass' MakeMove
  MultiStageClass::MakeMove();
}


void
IonDisplaceMoveClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar ("Sigma", Sigma));
  string ionSpeciesStr;
  assert (in.ReadVar("IonSpecies",   ionSpeciesStr));
  IonSpeciesNum  = Path.SpeciesNum(ionSpeciesStr);
  if (IonSpeciesNum == -1) {
    cerr << "Unrecogonized species, """ << ionSpeciesStr 
	 << """ in IonDisplaceMoveClass::Read.\n";
    abort();
  }
  assert (in.ReadVar ("NumToMove", NumParticlesToMove));
  in.ReadVar("WarpElectrons", WarpElectrons, false);

  IonDisplaceStage.Sigma = Sigma;
  IonDisplaceStage.IonSpeciesNum = IonSpeciesNum;
  IonDisplaceStage.WarpElectrons = WarpElectrons;

  // Construct action list
  IonDisplaceStage.Actions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange) 
    if (PathData.Actions.UseRPA)
      IonDisplaceStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    else
      IonDisplaceStage.Actions.push_back(&PathData.Actions.LongRange);

  if ((PathData.Actions.NodalActions(IonSpeciesNum)!=NULL)) {
    cerr << "IonDisplaceMove adding fermion node action for species " 
	 << ionSpeciesStr << endl;
    IonDisplaceStage.Actions.push_back
      (PathData.Actions.NodalActions(IonSpeciesNum));
  }
    
  // Now construct stage list
  Stages.push_back(&IonDisplaceStage);

  ActiveParticles.resize(NumParticlesToMove);
}
