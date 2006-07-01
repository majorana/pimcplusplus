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
    rhat.resize(N);
  }
  dVec zero(0.0);
  DeltaRions = zero;
  /// Now, choose a random displacement 
  for (int ptclIndex=0; ptclIndex<IonsToMove.size(); ptclIndex++) {
    int ptcl = IonsToMove(ptclIndex);
    Vec3 delta;
    PathData.Path.Random.CommonGaussianVec (Sigma, delta);
    DeltaRions(ptcl-first) = delta;

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

inline double det (const TinyMatrix<double,3,3> &A)
{
  return 
    A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1)) -
    A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0)) +
    A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
}



double
IonDisplaceStageClass::DoElectronWarp()
{
  cerr << "DoElectronWarp...\n";
  SpeciesClass &ionSpecies = Path.Species(IonSpeciesNum);
  int ionFirst = ionSpecies.FirstPtcl;
  int ionLast  = ionSpecies.LastPtcl;
  int N     = ionLast - ionFirst + 1;
  dVec disp;
  double dist;
  double logJProdForw = 0.0;
  double logJProdRev  = 0.0;
  // The jacobian matrices for the forward and reverse move
  TinyMatrix<double,3,3> J;
  Array<double,1> g(N);
  Array<Vec3,1> wgr(N);
  for (int si=0; si<Path.NumSpecies(); si++) {
    SpeciesClass &species = Path.Species(si);    
    if (species.lambda > 1.0e-10) {
      for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
	for (int elec=species.FirstPtcl; elec<=species.LastPtcl; elec++) {
	  SetMode (OLDMODE);
	  Weights = 0.0;
	  double totalWeight = 0.0;
	  for (int ion=ionFirst; ion <= ionLast; ion++) {
	    Path.DistDisp(slice, ion, elec, dist, disp);
	    double d4 = dist*dist*dist*dist;
	    Weights(ion-ionFirst) = 1.0/d4;
	    g(ion-ionFirst)      = -4.0/dist;
	    rhat(ion-ionFirst )   = (1.0/dist)*disp;
	    totalWeight += Weights(ion-ionFirst);
	  }
	  Weights *= (1.0/totalWeight);
	  SetMode (NEWMODE);
	  for (int ion=ionFirst; ion <= ionLast; ion++) 
	    Path(slice, elec) = Path(slice,elec) + 
	      Weights(ion-ionFirst)*DeltaRions(ion-ionFirst);

	  // Calculate Jacobian for forward move
	  J = 1.0, 0.0, 0.0,
  	      0.0, 1.0, 0.0,
	      0.0, 0.0, 1.0;
	  Vec3 sumwgr(0.0, 0.0, 0.0);
	  for (int j=0; j<N; j++) {
	    wgr(j) = Weights(j)*g(j)*rhat(j);
	    sumwgr += wgr(j);
	  }
	  for (int i=0; i<N; i++) {
	    // Vec3 dwdr = Weights(i)*(g(i)*rhat(i) - sumwgr);
	    Vec3 dwdr = (wgr(i) - Weights(i)*sumwgr);
	    for (int alpha=0; alpha<3; alpha++)
	      for (int beta=0; beta<3; beta++) 
		J(alpha,beta) += DeltaRions(i)[alpha] * dwdr[beta];
	  }
	  logJProdForw += log(det(J));
	}
      }
    }
  }
  // Repeat the loop for the reverse move to calculate the reverse
  // Jacobian 
  DeltaRions = -1.0*DeltaRions;
  for (int si=0; si<Path.NumSpecies(); si++) {
    SpeciesClass &species = Path.Species(si);    
    if (species.lambda > 1.0e-10) {
      for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
	for (int elec=species.FirstPtcl; elec<=species.LastPtcl; elec++) {
	  Weights = 0.0;
	  double totalWeight = 0.0;
	  for (int ion=ionFirst; ion <= ionLast; ion++) {
	    Path.DistDisp(slice, ion, elec, dist, disp);
	    double d4 = dist*dist*dist*dist;
	    Weights(ion-ionFirst) = 1.0/d4;
	    g(ion-ionFirst)      = -4.0/dist;
	    rhat(ion-ionFirst )   = (1.0/dist)*disp;
	    totalWeight += Weights(ion-ionFirst);
	  }
	  Weights *= (1.0/totalWeight);
	  
	  // Calculate Jacobian for forward move
	  J = 1.0, 0.0, 0.0,
  	      0.0, 1.0, 0.0,
	      0.0, 0.0, 1.0;
	  Vec3 sumwgr(0.0, 0.0, 0.0);
	  for (int j=0; j<N; j++) {
	    wgr(j) = Weights(j)*g(j)*rhat(j);
	    sumwgr += wgr(j);
	  }
	  for (int i=0; i<N; i++) {
	    Vec3 dwdr = (wgr(i) - Weights(i)*sumwgr);
	    for (int alpha=0; alpha<3; alpha++)
	      for (int beta=0; beta<3; beta++) 
		J(alpha,beta) += DeltaRions(i)[alpha] * dwdr[beta];
	  }
	  logJProdRev += log(det(J));
	}
      }
    }
  }
  cerr << "JProdForw = " << exp(logJProdForw) << endl;
  cerr << "JProdRev =  " << exp(logJProdRev)  << endl;
  cerr << "JProdRev/JProdForw = " << exp(logJProdRev-logJProdForw) << endl;
  //  return exp(logJProdRev-logJProdForw);
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
  for (int i=0; i<NumIonsToMove; i++) {
    int index = PathData.Path.Random.CommonInt(numLeft);
    vector<int>::iterator iter = ptclList.begin();
    IonDisplaceStage.IonsToMove(i) = ptclList[index];
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
  assert (in.ReadVar ("NumToMove", NumIonsToMove));
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

  IonDisplaceStage.IonsToMove.resize(NumIonsToMove);
  if (WarpElectrons) {
    //    IonDisplaceStage.Actions.push_back(&PathData.Actions.Kinetic);
    ActiveParticles.resize(Path.NumParticles());
    for (int i=0; i<Path.NumParticles(); i++)
      ActiveParticles(i) = i;
  }
  else {
    SpeciesClass &species = Path.Species(IonSpeciesNum);
    ActiveParticles.resize(species.NumParticles);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
      ActiveParticles(ptcl-species.FirstPtcl) = ptcl;
  }
}
