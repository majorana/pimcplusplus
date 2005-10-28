#include "CorrelatedBisectionBlock.h"
#include "StructureReject.h"

void CorrelatedBisectionBlockClass::Read(IOSectionClass &in)
{

  string permuteType, speciesName;
  //  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  assert (in.ReadVar ("name", Name));


  BisectionStages.resize(NumLevels);
  for (int level=0;level<NumLevels;level++)
    BisectionStages(NumLevels-level-1)=new BisectionStageClass(PathData,level,OutSection);


  
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  if    (PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION){
    HaveRefslice=true;
  }
 HaveRefslice = 
    ((PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION) &&
     (PathData.Actions.NodalActions(SpeciesNum) != NULL) &&
     (!PathData.Actions.NodalActions(SpeciesNum)->IsGroundState()));
  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    PermuteStage = new TablePermuteStageClass(PathData, SpeciesNum, NumLevels,
					      OutSection);
  else if (permuteType=="COUPLE")
    PermuteStage= new CoupledPermuteStageClass(PathData,SpeciesNum,NumLevels,
					       OutSection);
  else if (permuteType == "NONE") 
    PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels,
					   OutSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);

  LocalActions.push_back(&PathData.Actions.Kinetic);
  CommonActions.push_back(&PathData.Actions.Kinetic);
  LocalActions.push_back(&PathData.Actions.ShortRange);
  CommonActions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange){
    if (PathData.Actions.UseRPA){
      LocalActions.push_back(&PathData.Actions.LongRangeRPA);
      CommonActions.push_back(&PathData.Actions.LongRangeRPA);
    }
    else{
      LocalActions.push_back(&PathData.Actions.LongRange);
      CommonActions.push_back(&PathData.Actions.LongRange);
    }
  }
  for (int species=0;species<PathData.NumSpecies();species++)
    if ((PathData.Actions.NodalActions(species)!=NULL)) {
      cerr << "Adding fermion node action for species " 
	   << speciesName << endl;
      CommonActions.push_back(PathData.Actions.NodalActions(species));
    }
}


void CorrelatedBisectionBlockClass::ChooseTimeSlices()
{
  //  if (PathData.Path.Communicator.MyProc()==0)
    //    cerr<<"Choosing time slices"<<endl;
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  // do something special to avoid moving reference slice
  if (HaveRefslice &&
      Path.SliceOwner(Path.GetRefSlice()) == myProc) {
    cerr << "Avoid reference slice.\n";
    int bSlices = 1<<NumLevels;
    int myStart, myEnd;
    Path.SliceRange (myProc, myStart, myEnd);
    // refSlice is relative to my first slice
    int refSlice = Path.GetRefSlice() - myStart;
    int numSlices = Path.NumTimeSlices();
    assert(bSlices*2<numSlices);
    if (refSlice < bSlices) {  
      int numStarts = numSlices - bSlices - refSlice;
      Slice1 = refSlice + Path.Random.LocalInt(numStarts);
      Slice2 = Slice1 + bSlices;
      if (Slice2>=numSlices){
	cerr << "(Slice 2, numSlices): (" << Slice2 << ", " 
	     << numSlices << ")" << " "<< bSlices << " "
	     << numStarts << " " << refSlice << endl;
      }
      assert (Slice2 < numSlices);
    }
    else if (refSlice > (numSlices-1-bSlices) && refSlice<numSlices) {
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
  //  else if (PathData.Path.OpenPtcl) {
    

  //  }
  else {
    int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    Slice1 = PathData.Path.Random.LocalInt (numLeft);
    Slice2 = Slice1+sliceSep;
  }
  //  if (PathData.Path.Communicator.MyProc()==0)
    //    cerr<<"Ending Choosing time slices"<<endl;
}

void CorrelatedBisectionBlockClass::MakeMove()
{
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);
  PermuteStage->InitBlock();
  ActiveParticles.resize(1);
  Array<int,1> allParticles(PathData.Path.NumParticles());
  for (int i=0; i<allParticles.size(); i++)
    allParticles(i) = i;
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0) = -1;
    double logSampleProb=0.0;
    ///This is where we call the permtuations
    double prevActionChange=0.0;
    PermuteStage->Attempt(Slice1,Slice2,ActiveParticles,prevActionChange);
    for (int bIndex=0;bIndex<BisectionStages.size();bIndex++)
      logSampleProb += 
	log(BisectionStages(bIndex)->Sample(Slice1,Slice2,ActiveParticles));
    
    list<ActionBaseClass*>::iterator iter;
    PathData.Path.SetIonConfig(0);
    double localSAOld, localSANew, localSBOld, localSBNew;
    localSAOld = 0.0; localSANew = 0.0; localSBOld = 0.0; localSBNew =0.0;
    for (iter=LocalActions.begin();iter!=LocalActions.end();iter++){
      SetMode(NEWMODE);
      localSANew += 
	(*iter)->Action(Slice1,Slice2,ActiveParticles,0);
      SetMode(OLDMODE);
      localSAOld += 
	(*iter)->Action(Slice1,Slice2,ActiveParticles,0);
    }
    PathData.Path.SetIonConfig(1);
    for (iter=LocalActions.begin();iter!=LocalActions.end();iter++){
      SetMode(NEWMODE);
      localSBNew += 
	(*iter)->Action(Slice1,Slice2,ActiveParticles,0);
      SetMode(OLDMODE);
      localSBOld += 
	(*iter)->Action(Slice1,Slice2,ActiveParticles,0);
    }
    double localSBarOld=0.5*(localSAOld+localSBOld);
    double localSBarNew=0.5*(localSANew+localSBNew);
    double logAcceptProb= logSampleProb - localSBarNew + localSBarOld;
    bool toAcceptLocal = logAcceptProb>=log(PathData.Path.Random.Local());
    
    ///We compute here over all time slices becasue we need the total
    ///action.  Need to fix this to update the action dynamically. 
    int s1 = 0;
    int s2 = PathData.Path.NumTimeSlices()-1;

    /// Calculate new and old total sA and SB
    double totalSAOld, totalSANew, totalSBOld, totalSBNew;
    totalSAOld = 0.0; totalSANew = 0.0; totalSBOld = 0.0; totalSBNew = 0.0;

    //// MUST MOVE JOIN OUT OF THE WAY FOR THIS TO WORK
    SetMode(OLDMODE);
    for (iter=CommonActions.begin(); iter!=CommonActions.end(); iter++) 
      totalSBOld+=(*iter)->Action(s1, s2, allParticles, 0);
    SetMode(NEWMODE);    
    for (iter=CommonActions.begin(); iter!=CommonActions.end(); iter++) 
      totalSBNew+=(*iter)->Action(s1, s2, allParticles, 0);

    PathData.Path.SetIonConfig(0);
    SetMode(OLDMODE);
    for (iter=CommonActions.begin(); iter!=CommonActions.end(); iter++) 
      totalSAOld+=(*iter)->Action(s1, s2, allParticles, 0);
    SetMode(NEWMODE);    
      for (iter=CommonActions.begin(); iter!=CommonActions.end(); iter++) 
	totalSANew+=(*iter)->Action(s1, s2, allParticles, 0);

    double deltaSOld = 0.5*(totalSAOld - totalSBOld);
    double deltaSNew = 0.5*(totalSANew - totalSBNew);

    /// Only add my contribution to the sum if I've accepted locally
    if (toAcceptLocal) {
      deltaSNew  = PathData.Path.Communicator.AllSum(deltaSNew);
      totalSANew = PathData.Path.Communicator.AllSum(totalSANew);
      totalSBNew = PathData.Path.Communicator.AllSum(totalSBNew);
    }
    else { /// otherwise just add the old action 
      deltaSNew  = PathData.Path.Communicator.AllSum(deltaSOld);
      totalSANew = PathData.Path.Communicator.AllSum(totalSAOld);
      totalSBNew = PathData.Path.Communicator.AllSum(totalSBOld);

    }

    deltaSOld  = PathData.Path.Communicator.AllSum(deltaSOld);
    totalSAOld = PathData.Path.Communicator.AllSum(totalSAOld);
    totalSBOld = PathData.Path.Communicator.AllSum(totalSBOld);
    double commonAcceptProb = cosh(deltaSNew)/cosh(deltaSOld);
    
//     double commonAcceptProb = exp(logSampleProb) * 
//       (exp(-totalSANew)+exp(-totalSBNew))/(exp(-totalSAOld)+exp(-totalSBOld));

    bool toAcceptCommon = 
      log(commonAcceptProb)>=log(PathData.Path.Random.Common());
    
    double wA, wB;
    if (toAcceptCommon) {
      /// Accept
//       wA = exp(-totalSANew+0.5*(totalSANew+totalSBNew))/
// 	(2.0*cosh(deltaSNew));
      wA = exp(-deltaSNew)/(2.0*cosh(deltaSNew));
      //       wB = exp(-totalSBNew+0.5*(totalSANew+totalSBNew))/
      // 	(2.0*cosh(deltaSNew));
      wB = exp(+deltaSNew)/(2.0*cosh(deltaSNew));
      //      wB = 1.0-wA;
      assert (fabs(wA+wB-1.0) < 1.0e-12);
      if (toAcceptLocal)
	Accept();
      else
	Reject();
    }
    else {
      /// Reject
//       wA = exp(-totalSAOld+0.5*(totalSAOld+totalSBOld))/
// 	(2.0*cosh(deltaSOld));
//       wB = exp(-totalSBOld+0.5*(totalSAOld+totalSBOld))/
// 	(2.0*cosh(deltaSOld));
      wA = exp(-deltaSOld)/(2.0*cosh(deltaSOld));
      wB = exp(+deltaSOld)/(2.0*cosh(deltaSOld));
      assert (fabs(wA+wB-1.0) < 1.0e-12);
      wB = 1.0-wA;
      Reject();
    }

    double energyA, energyB, totalSA, totalSB;
    energyA = energyB = 0.0;
    totalSA = totalSB = 0.0;

    PathData.Path.SetIonConfig(0);
    for (iter=CommonActions.begin(); iter!=CommonActions.end(); iter++) {
      energyA += (*iter)->d_dBeta(s1,s2,0);
      totalSA += (*iter)->Action(s1, s2, allParticles, 0);
    }

    PathData.Path.SetIonConfig(1);
    for (iter=CommonActions.begin(); iter!=CommonActions.end(); iter++) {
      energyB += (*iter)->d_dBeta(s1,s2,0);
      totalSB += (*iter)->Action(s1, s2, allParticles,0);
    }

    energyA = PathData.Path.Communicator.AllSum(energyA);
    energyB = PathData.Path.Communicator.AllSum(energyB);
    totalSA = PathData.Path.Communicator.AllSum(totalSA);
    totalSB = PathData.Path.Communicator.AllSum(totalSB);
    energyA /= PathData.Path.TotalNumSlices;
    energyB /= PathData.Path.TotalNumSlices;

    //    cerr << "wA = " << wA << "    wB = " << wB << endl;
    PathData.Path.ActionA.push_back(totalSA);
    PathData.Path.WeightA.push_back(wA);
    PathData.Path.EnergyA.push_back(energyA);
    
    PathData.Path.ActionB.push_back(totalSB);
    PathData.Path.WeightB.push_back(wB);
    PathData.Path.EnergyB.push_back(energyB);

    StepNum++;
  }
  double wASASum = 0.0;
  double wBSBSum = 0.0;
  double wASum   = 0.0;
  double wBSum   = 0.0;
  double wAEASum   = 0.0;
  double wBEBSum   = 0.0;
  if ((StepNum % 100) == 99) {
    for (int i=0; i<PathData.Path.ActionA.size(); i++) {
      wASASum += (PathData.Path.ActionA[i]*PathData.Path.WeightA[i]);
      wAEASum += (PathData.Path.EnergyA[i]*PathData.Path.WeightA[i]);
      wASum   += PathData.Path.WeightA[i];
      wBSBSum += (PathData.Path.ActionB[i]*PathData.Path.WeightB[i]);
      wBEBSum += (PathData.Path.EnergyB[i]*PathData.Path.WeightB[i]);
      wBSum   += PathData.Path.WeightB[i];
      
    }
    perr << "SA = " << (wASASum/wASum) << "  " 
	 << "SB = " << (wBSBSum/wBSum) 
	 << "  deltaS = " << (wASASum/wASum - wBSBSum/wBSum) << endl;
    perr << "delta E = " << (wAEASum/wASum - wBEBSum/wBSum) << endl;
    perr << "AcceptRatio = " << ((double)NumAccepted/NumMoves) << endl;
//     fprintf (EAout, "%1.12e\n", (wAEASum/wASum));
//     fprintf (EBout, "%1.12e\n", (wBEBSum/wBSum));
//     fprintf (SAout, "%1.12e\n", (wASASum/wASum));
//     fprintf (SBout, "%1.12e\n", (wBSBSum/wBSum));
//     fflush (EAout);
//     fflush (EBout);
//     fflush (SAout);
//     fflush (SBout);

    wAEAvar.Write(wAEASum);
    wBEBvar.Write(wBEBSum);
    wASAvar.Write(wASASum);
    wBSBvar.Write(wBSBSum);
    wAvar.Write(wASum);
    wBvar.Write(wBSum);
    wAvar.Flush();

    PathData.Path.ActionA.clear();
    PathData.Path.WeightA.clear();
    PathData.Path.EnergyA.clear();
    PathData.Path.ActionB.clear();
    PathData.Path.WeightB.clear();
    PathData.Path.EnergyB.clear();

  }
}
