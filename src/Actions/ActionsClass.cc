/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "ActionsClass.h"
#include "../PathDataClass.h"
#include <Common/Ewald/OptimizedBreakup.h>
#include <Common/Integration/GKIntegration.h>
#include <Common/IO/FileExpand.h>


///Actionsclass. Stores all the actsion
void 
ActionsClass::Read(IOSectionClass &in)
{ 
  cerr<<"Starting Actions Read"<<endl;
#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  assert(in.ReadVar ("tau", Path.tau));
  assert(in.ReadVar ("MaxLevels", MaxLevels));
  assert(in.ReadVar ("NumImages", NumImages));
  Kinetic.SetNumImages (NumImages);
  KineticSphere.SetNumImages(NumImages);
  Mu.Read(in);
	bool doMolRead = false;
	if(in.ReadVar("InitMoleculeInteractions",doMolRead)){
		MoleculeInteractions.Read(in);
		MoleculeInteractions.SetNumImages(NumImages);
	}
  bool doLRCRead = false;
	if(in.ReadVar("InitLongRangeCoulomb",doLRCRead)){
		LongRangeCoulomb.Read(in);
	}
	bool doQBoxRead = false;
	if(in.ReadVar("InitQBoxAction",doQBoxRead)){
		QBoxAction.Read(in);
	}

#ifdef ORDER_N_FERMIONS
  VariationalPI.Read(in);
  TruncatedInverse.Read(in);
#endif
  //  VariationalPI.Read(in);
  perr << "MaxLevels = " << MaxLevels << endl;
  bool checkJosephson=false;
  in.ReadVar("Josephson",checkJosephson);
  if (checkJosephson){
    Josephson.Read(in);
    Hermele.Read(in);
    DualHermele.Read(in);
  }
  bool usePairAction=false;
  in.ReadVar("Paired",usePairAction);
  if (usePairAction){
    PairFixedPhase.Read(in);
  }
  
  if (!in.ReadVar ("UseRPA", UseRPA))
    UseRPA = false;
  if (UseRPA) 
    perr << "Using RPA for long range action.\n";
  else      
    perr << "Not using RPA for long range action.\n";

cerr << " going to read pair actions" << endl;
  Array<string,1> PAFiles;
  assert (in.ReadVar ("PairActionFiles", PAFiles));
  Array<string,1> SpecificHeatPAFiles;
  bool readSpecificHeatFiles=false;
  if (in.ReadVar("SpecificHeatPAFiles",SpecificHeatPAFiles))
    readSpecificHeatFiles=true;
  

  int numPairActions = PAFiles.size();
cerr << "Looking for " << numPairActions << endl;
  PairArray.resize(numPairActions);
  if (readSpecificHeatFiles)
    SpecificHeatPairArray.resize(SpecificHeatPAFiles.size());
  PairMatrix.resize(Path.NumSpecies(),Path.NumSpecies());
  // Initialize to a nonsense value so we can later check in the table
  // element was filled in.
  for (int i=0; i<Path.NumSpecies(); i++)
    for (int j=0; j<Path.NumSpecies(); j++)
      PairMatrix(i,j) = (PairActionFitClass*)NULL;
  // Read pair actions files
cerr << "declaring IOSectionClass...";
  IOSectionClass PAIO;
cerr << " done" << endl;

  for (int i=0; i<numPairActions; i++) {
    // Allow for tilde-expansion in these files
    string name = ExpandFileName(PAFiles(i));
    assert(PAIO.OpenFile (name));
cerr << i << ": reading " << name << endl;
    PairArray(i) = ReadPAFit (PAIO, Path.tau, MaxLevels);
    bool paUsed=false;
    for (int spec1=0;spec1<Path.NumSpecies();spec1++)
      for (int spec2=spec1;spec2<Path.NumSpecies();spec2++) 
	if (((Path.Species(spec1).Type==PairArray(i)->Particle1.Name)&&
	     (Path.Species(spec2).Type==PairArray(i)->Particle2.Name)) ||
	    ((Path.Species(spec2).Type==PairArray(i)->Particle1.Name)&&
	     (Path.Species(spec1).Type==PairArray(i)->Particle2.Name))) {
	  if (PairMatrix(spec1,spec2) != NULL) {
	    perr << "More than one pair action for species types (" 
		 << PairArray(i)->Particle1.Name << ", "
		 << PairArray(i)->Particle2.Name << ")." << endl;
	    exit(-1);
	  }
	  perr << "Found PAfile for pair (" 
	       << Path.Species(spec1).Name << ", "
	       << Path.Species(spec2).Name << ")\n";
	  PairMatrix(spec1,spec2) = PairArray(i);
	  PairMatrix(spec2,spec1) = PairArray(i);
	  paUsed = true;
	}
    if (!paUsed) {
      perr << "Warning:  Pair action for species types (" 
	   << PairArray(i)->Particle1.Name << ", "
 	   << PairArray(i)->Particle1.Name << ") not used.\n";
    }

    PAIO.CloseFile();
  }

  if (readSpecificHeatFiles){
    cerr<<"I READ SPECIFIC HEAT FILES"<<endl;
    assert((in.ReadVar("TauValues",TauValues)));
    assert(TauValues.size()==SpecificHeatPAFiles.size());
    for (int i=0; i<SpecificHeatPAFiles.size() ; i++) {
      // Allow for tilde-expansion in these files
      string name = ExpandFileName(SpecificHeatPAFiles(i));
      assert(PAIO.OpenFile (name));
      SpecificHeatPairArray(i) = ReadPAFit (PAIO, TauValues(i), MaxLevels);	
      PAIO.CloseFile();
    }
  }

   
  // Now check to make sure all PairActions that we need are defined.
  for (int species1=0; species1<Path.NumSpecies(); species1++)
    for (int species2=0; species2<Path.NumSpecies(); species2++)
      if (PairMatrix(species1,species2) == NULL) {
	if ((species1 != species2) || 
	    (Path.Species(species1).NumParticles > 1)) {
	  perr << "We're missing a PairAction for species1 = "
	       << Path.Species(species1).Name << " and species2 = "
	       << Path.Species(species2).Name << endl;
	  exit(1);
	}
      }
  
  in.ReadVar("UseLongRange", UseLongRange);
  if (HaveLongRange()) {
    perr << "*** Using long-range/short-range breakup. ***\n";
    assert (in.ReadVar("UseBackground", LongRange.UseBackground));
    LongRangePot.UseBackground = LongRange.UseBackground;
    LongRangeRPA.UseBackground = LongRange.UseBackground;
  }

//   if (longRange){
//     LongRange.Init(in);
//     if (UseRPA)
//       LongRangeRPA.Init(in);
//   }

//  Create nodal action objects
//   NodalActions.resize(PathData.Path.NumSpecies());
//   for (int spIndex=0; spIndex<PathData.Path.NumSpecies(); spIndex++) {
//     SpeciesClass &species = PathData.Path.Species(spIndex);
//     if (species.GetParticleType() == FERMION) {
//       if (species.NodeType == "FREE") 
// 	NodalActions (spIndex) = new FreeNodalActionClass(PathData, spIndex);
//       else {
// 	cerr << "Unrecognized node type " << species.NodeType << ".\n";
// 	exit(EXIT_FAILURE);
//       }
//       NodalActions(spIndex)->Read(in);
//     }
//     else
//       NodalActions(spIndex) = NULL;
//   }
  
  ReadNodalActions (in);

//   // Now create nodal actions for Fermions
//   NodalActions.resize(PathData.Path.NumSpecies());
//   for (int species=0; species<PathData.Path.NumSpecies(); species++) 
//     if (PathData.Path.Species(species).GetParticleType() == FERMION)
//       NodalActions(species) = new FPNodalActionClass(PathData, species);
//     else
//       NodalActions(species) = NULL;
  
  perr << "Finished reading the action.\n"; 

  if (in.OpenSection("StructureReject")){
    StructureReject.Read(in);
    in.CloseSection();
  }
  if (in.OpenSection("Tether")){
    Tether.Read(in);
    in.CloseSection();
  }
  ///Reading in information for David long range action
  if (PathData.Path.DavidLongRange){
    DavidLongRange.Read(in);
  }
  if (PathData.Path.OpenPaths)
    OpenLoopImportance.Read(in);
}

/// Read in the nodal actions.
/// This should only be called after the PairActions have been read.
void
ActionsClass::ReadNodalActions(IOSectionClass &in)
{
  int numNodeSections=in.CountSections("NodalAction");
  NodalActions.resize (PathData.Path.NumSpecies());
  NodalActions = NULL;
  for (int nodeSection=0; nodeSection<numNodeSections; nodeSection++) {
    in.OpenSection("NodalAction", nodeSection);
    string type, speciesString;
    assert (in.ReadVar ("Type", type));
    if (type == "FREE") {
      assert (in.ReadVar("Species", speciesString));
      int species = PathData.Path.SpeciesNum(speciesString);

      FreeNodalActionClass &nodeAction = 
	*(new FreeNodalActionClass (PathData, species));
      nodeAction.Read(in);
      NodalActions(species) = &nodeAction;
    }
    else if (type == "GROUNDSTATE") {
      GroundStateClass &groundState = *new GroundStateClass(PathData);
      groundState.Read (in);
      NodalActions(groundState.UpSpeciesNum) = 
	new GroundStateNodalActionClass 
	(PathData, groundState, groundState.UpSpeciesNum);
      NodalActions(groundState.DownSpeciesNum) = 
	new GroundStateNodalActionClass 
	(PathData, groundState, groundState.DownSpeciesNum);
      NodalActions(groundState.IonSpeciesNum) = 
	new GroundStateNodalActionClass 
	(PathData, groundState, groundState.IonSpeciesNum);
    }
     else if (type == "VARIATIONAL") {
	 NodalActions.resizeAndPreserve(2);
       NodalActions(0) = 
 	&VariationalPI;
     }
    else if (type == "TRUNCATED") {
	 NodalActions.resizeAndPreserve(2);
      NodalActions(1) = 
	&TruncatedInverse;
    }
    else if (type == "FIXEDPHASE") {
      FixedPhaseClass &fixedPhaseA = *new FixedPhaseClass(PathData);
      fixedPhaseA.Read (in);
      FixedPhaseClass &fixedPhaseB = 
	PathData.Path.UseCorrelatedSampling() ? *new FixedPhaseClass(PathData) : fixedPhaseA;
      if (PathData.Path.UseCorrelatedSampling())
	fixedPhaseB.Read (in);
      NodalActions(fixedPhaseA.UpSpeciesNum) = 
	new FixedPhaseActionClass 
	(PathData, fixedPhaseA, fixedPhaseB, fixedPhaseA.UpSpeciesNum);
      NodalActions(fixedPhaseA.DownSpeciesNum) = 
	new FixedPhaseActionClass 
	(PathData, fixedPhaseA, fixedPhaseB, fixedPhaseA.DownSpeciesNum);
      NodalActions(fixedPhaseA.IonSpeciesNum) = 
	new FixedPhaseActionClass 
	(PathData, fixedPhaseA, fixedPhaseB, fixedPhaseA.IonSpeciesNum);
    }
    in.CloseSection();
  }
  cerr<<"Ending actions read"<<endl;
}



void
ActionsClass::Energy (double& kinetic, double &dUShort, double &dULong, 
		      double &node, double &vShort, double &vLong)
//void ActionsClass::Energy(map<double>& Energies)
{
	//double kinetic, dUShort, dULong, node, vShort, vLong;
  bool doLongRange = HaveLongRange() && UseLongRange;
  int M = PathData.Path.NumTimeSlices()-1;
  kinetic = Kinetic.d_dBeta (0, M, 0);
  if (PathData.Path.OrderN)
    dUShort=ShortRangeOn.d_dBeta(0,M,0);
  else
    dUShort = ShortRange.d_dBeta (0, M, 0);
  dULong=0.0;
  if (doLongRange){
    if (UseRPA)
      dULong = LongRangeRPA.d_dBeta (0, M, 0);
    else
      dULong = LongRange.d_dBeta (0, M, 0);
  }
  if (PathData.Path.DavidLongRange)
    dULong = DavidLongRange.d_dBeta(0,M,0);
  node = 0.0;
  for (int species=0; species<PathData.Path.NumSpecies(); species++)
    if (NodalActions(species) != NULL)
      node += NodalActions(species)->d_dBeta(0, M, 0);

  vShort  = 0.0; vLong   = 0.0;  
  for (int slice=0; slice <= M; slice++) {
    double factor = ((slice==0)||(slice==M)) ? 0.5 : 1.0;
    vShort += factor * ShortRangePot.V(slice);
    if (doLongRange)
      vLong  += factor *  LongRangePot.V(slice);
  }
	//Energies["kinetic"] = kinetic;
	//Energies["dUShort"] = dUShort;
	//Energies["dULong"] = dULong;
	//Energies["node"] = node;
	//Energies["vShort"] = vShort;
	//Energies["vLong"] = vLong;
}



void
ActionsClass::GetActions (double& kinetic, double &UShort, double &ULong, 
			  double &node)
{
  bool doLongRange = HaveLongRange() && UseLongRange;
  Array<int,1> activePtcls(PathData.Path.NumParticles());
  for (int i=0; i<PathData.Path.NumParticles(); i++)
    activePtcls(i) = i;
  
  int M = PathData.Path.NumTimeSlices()-1;
  kinetic = Kinetic.Action (0, M, activePtcls, 0);
  UShort = ShortRange.Action (0, M, activePtcls, 0);
  ULong=0.0;
  if (doLongRange){
    if (UseRPA)
      ULong = LongRangeRPA.Action (0, M, activePtcls, 0);
    else if (PathData.Path.DavidLongRange)
      ULong = DavidLongRange.Action(0,M, activePtcls, 0);
    else
      ULong = LongRange.Action (0, M, activePtcls, 0);
  }
  node = 0.0;
  for (int species=0; species<PathData.Path.NumSpecies(); species++)
    if (NodalActions(species) != NULL)
      node += NodalActions(species)->Action(0, M, activePtcls, 0);
}


//   PathClass &Path = PathData.Path;
//   for (int link=0; link<Path.NumTimeSlices()-1; link++) {    
//     for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
//       int specNum1 = Path.ParticleSpeciesNum(ptcl1);
//       SpeciesClass &spec1 = Path.Species(specNum1);
//       if (spec1.lambda != 0.0) {
// 	// Compute kinetic energy
// 	/// Add constant part to kinetic part of energy
// 	kinetic += (NDIM*0.5)/tau;
// 	// Now do spring part
// 	double fourLambdaTauInv = 1.0/(4.0*spec1.lambda*tau);
// 	dVec vel;
// 	vel = Path.Velocity (link, link+1, ptcl1);
// 	double Z = 1.0;
// 	dVec gaussSum = 0.0;
// 	dVec numSum = 0.0;
// 	for (int dim=0; dim<NDIM; dim++) {
// 	  for (int image=-NumImages; image<=NumImages; image++) {
// 	    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
// 	    double dist2OverFLT = dist*dist*fourLambdaTauInv;
// 	    double expPart = exp(-dist2OverFLT);
// 	    gaussSum[dim] += expPart;
// 	    numSum[dim]   += dist2OverFLT*expPart/tau;
// 	  }
// 	  Z *= gaussSum[dim];
// 	}
      
// 	double scalarnumSum = 0.0;
// 	for (int dim=0;dim<NDIM;dim++){
// 	  dVec numProd=1.0;
// 	  for (int dim2=0;dim2<NDIM;dim2++)
// 	    if (dim2!=dim)
// 	      numProd[dim] *= gaussSum[dim2];
// 	    else 
// 	      numProd[dim] *=  numSum[dim2];
// 	  scalarnumSum += numProd[dim];
// 	}
// 	kinetic += scalarnumSum/Z; 
//       }
    
//       // Now do short-range part of energy
//       for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
// 	int specNum2 = Path.ParticleSpeciesNum(ptcl2);
// 	dVec r, rp;
// 	double rmag, rpmag;
// 	Path.DistDisp(link, link+1, ptcl1, ptcl2, rmag,rpmag,r,rp); 
	
// 	double s2 = dot(r-rp, r-rp);
// 	double q = 0.5*(rmag+rpmag);
// 	double z = (rmag-rpmag);
	
// 	PairActionFitClass &pa = *PairMatrix(specNum1,specNum2);
// 	duShort += pa.dU(q, z, s2, 0);
// 	// Subtract off long-range part from short-range action
// 	if (pa.IsLongRange())
// 	  duShort -= 0.5*(pa.dUlong(0)(rmag)+pa.dUlong(0)(rpmag));
//       }
//     }
//   }

//   if (UseRPA)
//     dULong = LongRangeRPA.d_dBeta (0, Path.NumTimeSlices()-1, 0);
//   else
//     dULong = LongRange.d_dBeta (0, Path.NumTimeSlices()-1, 0);


//    // Now, calculate potential
//    for (int slice=0; slice < Path.NumTimeSlices; slice++) {
//      double factor = ((slice==0)||(slice==Path.NumTimesSlices()-1)) ? 0.5 : 1.0;
//      for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
//        int specNum1 = Path.ParticleSpeciesNum(ptcl1);
//        for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
// 	 int specNum2 = Path.ParticleSpeciesNum(ptcl2);
// 	 double dist;
// 	 dVec disp;
// 	 Path.DistDisp (slice, ptcl1, ptcl2, dist, disp);
// 	 PairActionFitClass &pa = *PairMatrix(specNum1,specNum2);
// 	 vShort += factor * pa.V(dist);
// 	 if (pa.IsLongRange())
// 	   vShort -= factor * pa.Vlong(dist);
//        }
//      }
//    }

     
	 

// }


Potential&
ActionsClass::GetPotential (int species1, int species2)
{
  return *(PairMatrix(species1, species2)->Pot);
}


void 
ActionsClass::ShiftData (int slicesToShift)
{
  OpenLoopImportance.ShiftData(slicesToShift);
  StructureReject.ShiftData(slicesToShift);
  ShortRange.ShiftData(slicesToShift);
  ShortRangeOn.ShiftData(slicesToShift);
  ShortRangeApproximate.ShiftData(slicesToShift);
  ShortRangePrimitive.ShiftData(slicesToShift);
  LongRange.ShiftData(slicesToShift);
  LongRangeRPA.ShiftData(slicesToShift);
  DavidLongRange.ShiftData(slicesToShift);
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i)!=NULL)
      NodalActions(i)->ShiftData(slicesToShift);
}


void 
ActionsClass::AcceptCopy (int startSlice, int endSlice,
			  const Array<int,1> &activeParticles)
{
#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  bool activeSpecies[Path.NumSpecies()];
  for (int i=0; i<Path.NumSpecies(); i++)
    activeSpecies[i] = false;
  for (int pi=0; pi<activeParticles.size(); pi++)
    activeSpecies[Path.ParticleSpeciesNum(activeParticles(pi))] = true;

  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL && activeSpecies[i])
      NodalActions(i)->AcceptCopy (startSlice, endSlice);

  QBoxAction.AcceptCopy(startSlice, endSlice);
}


void 
ActionsClass::RejectCopy (int startSlice, int endSlice,
			       const Array<int,1> &activeParticles)
{
#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  bool activeSpecies[Path.NumSpecies()];
  for (int i=0; i<Path.NumSpecies(); i++)
    activeSpecies[i] = false;
  for (int pi=0; pi<activeParticles.size(); pi++)
    activeSpecies[Path.ParticleSpeciesNum(activeParticles(pi))] = true;

  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL && activeSpecies[i])
      NodalActions(i)->RejectCopy (startSlice, endSlice);

  QBoxAction.RejectCopy(startSlice, endSlice);
}

void
ActionsClass::Init()
{
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->Init();
}


bool ActionsClass::HaveLongRange()
{
  bool longRange = false;
  for (int i=0; i<PairArray.size(); i++)
    longRange = longRange || PairArray(i)->IsLongRange();
  return (UseLongRange && longRange);
}

void 
ActionsClass::Setk(Vec3 k)
{
  for (int i=0; i<PathData.Path.NumSpecies(); i++) {
    if (NodalActions(i) != NULL)
      NodalActions(i)->Setk(k);
  }
}

void
ActionsClass::WriteInfo(IOSectionClass &out)
{
  /// If we have nodal actions, have them write any pertinent info to
  /// the output file.
  bool haveNodeActions = false;
  for (int i=0; i<PathData.Path.NumSpecies(); i++) 
    if (NodalActions(i) != NULL)
      haveNodeActions = true;
  if (haveNodeActions) {
    if (PathData.IntraComm.MyProc() == 0) {
      out.NewSection("NodalActions");
      for (int i=0; i<PathData.Path.NumSpecies(); i++) {
	out.NewSection ("NodeAction");
	if (NodalActions(i) == NULL) 
	  out.WriteVar("Type", "NONE");
	else 
	  NodalActions(i)->WriteInfo(out);
	
	out.CloseSection();
      }
      out.CloseSection();
    }
  }
}


void
ActionsClass::GetForces(const Array<int,1> &ptcls, 
			Array<dVec,1> &Fshort, Array<dVec,1> &Flong)
{
  //Move the join to the end so we don't have to worry about
  //permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  assert (Fshort.size() == ptcls.size());
  assert (Flong.size() == ptcls.size());
  Array<dVec,1> FtmpShort(Fshort.size());
  Array<dVec,1> FtmpLong(Flong.size());
  dVec zero; 
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  FtmpShort = zero;
  FtmpLong = zero;
  /// Calculate short and long-range pair actions for now.
  ShortRange.GradAction(0, Path.NumTimeSlices()-1, ptcls, 0, FtmpShort);
  LongRange.GradAction(0, Path.NumTimeSlices()-1, ptcls, 0, FtmpLong);
  /// Answer must be divided by beta.
  double beta = Path.TotalNumSlices * Path.tau;
  Fshort -= (1.0/beta)*FtmpShort;
  Flong -= (1.0/beta)*FtmpLong;
}


void
ActionsClass::GetForcesFD(const Array<int,1> &ptcls, Array<dVec,1> &F)
{
  const double eps = 1.0e-6;
#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  Array<dVec,1> Ftmp(F.size());
  dVec zero;
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  Ftmp = zero;
  Array<int,1> onePtcl(1);
  for (int pi=0; pi < ptcls.size(); pi++) {
    int ptcl = ptcls(pi);
    onePtcl(0) = ptcl;
    dVec savePos = Path(0, ptcl);
    dVec uPlus, uMinus;
    for (int dim=0; dim<NDIM; dim++) {
      for (int slice=0; slice<Path.NumTimeSlices(); slice++)
	Path(slice,ptcl)[dim] = savePos[dim] + eps;
      uPlus[dim] = ShortRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      uPlus[dim] += LongRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      for (int slice=0; slice<Path.NumTimeSlices(); slice++)
	Path(slice,ptcl)[dim] = savePos[dim] - eps;
      uMinus[dim] = ShortRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      uMinus[dim] += LongRange.Action(0, Path.NumTimeSlices()-1, onePtcl, 0);
      for (int slice=0; slice<Path.NumTimeSlices(); slice++)
	Path(slice,ptcl)[dim] = savePos[dim];
    }
    Ftmp(pi) = (0.5/eps)*(uPlus-uMinus);
  }
  double beta = Path.TotalNumSlices * Path.tau;
  F -= (1.0/beta)*Ftmp;
}


void
ActionsClass::UpdateNodalActions()
{
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->Update();
}

void
ActionsClass::MoveJoin (int oldJoinPos, int newJoinPos)
{
  // Currently, only some nodal actions actually need their MoveJoin called.
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->MoveJoin (oldJoinPos, newJoinPos);
}
