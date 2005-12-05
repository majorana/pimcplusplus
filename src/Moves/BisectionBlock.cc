#include "EmptyStage.h"
#include "BisectionBlock.h"
#include "StructureRejectStage.h"
#include "CouplingStage.h"
#include "WormPermuteStage.h"
#include "NoPermuteStage.h"


void BisectionBlockClass::Read(IOSectionClass &in)
{
  //  bool orderNBosons;
  string permuteType, speciesName;
  //  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  assert (in.ReadVar ("name", Name));
  //  in.ReadVar("OrderNBosons",orderNBosons);
  bool addStructureRejectStage=false;
  in.ReadVar("StructureReject",addStructureRejectStage);
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
					      IOSection);
  else if (permuteType=="COUPLE")
    PermuteStage= new CoupledPermuteStageClass(PathData,SpeciesNum,NumLevels,
					       IOSection);
  else if (permuteType == "NONE") 
    PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels,
					   IOSection);
  else if (permuteType=="OPEN")
    PermuteStage = new WormPermuteStage2Class(PathData,SpeciesNum,NumLevels,
					     IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);
  Stages.push_back (PermuteStage);
  if (permuteType=="OPEN"){
    EmptyStageClass *newStage=new EmptyStageClass(PathData,IOSection);
    newStage->Read(in);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    Stages.push_back(newStage);
  }
//   if (PathData.Path.Random.Local()>0.5)
//     PathData.Path.ExistsCoupling=1;
//   else
//     PathData.Path.ExistsCoupling=0;

  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level,
							     IOSection);
    newStage->TotalLevels=NumLevels;
    newStage->Actions.push_back(&PathData.Actions.Kinetic);

    
    if (PathData.Path.OpenPaths && level==0)
      newStage->Actions.push_back(&PathData.Actions.OpenLoopImportance);
    if (PathData.Path.OpenPaths && level>0)
      newStage->Actions.push_back(&PathData.Actions.ShortRangeApproximate);
    else if (PathData.Path.OrderN)
      newStage->Actions.push_back(&PathData.Actions.ShortRangeOn);
    else
      newStage->Actions.push_back(&PathData.Actions.ShortRange);
    if (level == 0) {
      if (addStructureRejectStage)
	newStage->Actions.push_back(&PathData.Actions.StructureReject);
      ///If it's David's long range class then do this
      if (PathData.Path.DavidLongRange){
	newStage->Actions.push_back(&PathData.Actions.DavidLongRange);
      }
      else if (PathData.Actions.HaveLongRange()){
	if (PathData.Actions.UseRPA)
	  newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
	else
	  newStage->Actions.push_back(&PathData.Actions.LongRange);
      }
      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	cerr << "Adding fermion node action for species " 
	     << speciesName << endl;
	newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
//       for (int i=0; i<PathData.Actions.NodalActions.size(); i++)
// 	newStage->Actions.push_back(PathData.Actions.NodalActions(i));
    }
    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }
  // Add the second stage of the permutation step
  Stages.push_back (PermuteStage);

//   ///HACK! Addding a stage that will reject the move if the structure
//   //factor gets too large
  bool useStructureRejectStage=false;
  in.ReadVar("StructureReject",useStructureRejectStage);
  if (useStructureRejectStage){
    StructureRejectStageClass* structureReject =
      new StructureRejectStageClass(PathData,in,IOSection);
    structureReject->Read(in);
    Stages.push_back(structureReject);
  }
  // Add the nodal action stage, if necessary
}


void BisectionBlockClass::ChooseTimeSlices()
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

void BisectionBlockClass::MakeMove()
{
  cerr<<"I am running bisection block"<<endl;
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);

  //  cerr<<"I am binning from "<<Slice1<<" to "<<Slice2<<"in bisection block"<<endl;
  if (PathData.Path.OrderN){
    for (int slice=Slice1;slice<=Slice2;slice++)
      PathData.Path.Cell.BinParticles(slice);
  }

  ((PermuteStageClass*)PermuteStage)->InitBlock();
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0)=-1;
    MultiStageClass::MakeMove();
  }
}
