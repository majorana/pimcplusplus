#include "ActionsClass.h"
#include "../PathDataClass.h"
#include <Common/Ewald/OptimizedBreakup.h>
#include <Common/Integration/GKIntegration.h>
#include <Common/IO/FileExpand.h>


///Actionsclass. Stores all the actsion
void 
ActionsClass::Read(IOSectionClass &in)
{ 
  PathClass &Path=PathData.Path;
  assert(in.ReadVar ("tau", Path.tau));
  assert(in.ReadVar ("MaxLevels", MaxLevels));
  assert(in.ReadVar ("NumImages", NumImages));
  Kinetic.SetNumImages (NumImages);
  KineticSphere.SetNumImages(NumImages);
  perr << "MaxLevels = " << MaxLevels << endl;

  if (!in.ReadVar ("UseRPA", UseRPA))
    UseRPA = false;
  if (UseRPA) 
    perr << "Using RPA for long range action.\n";
  else      
    perr << "Not using RPA for long range action.\n";


  Array<string,1> PAFiles;
  assert (in.ReadVar ("PairActionFiles", PAFiles));
  int numPairActions = PAFiles.size();
  PairArray.resize(numPairActions);
  PairMatrix.resize(Path.NumSpecies(),Path.NumSpecies());
  // Initialize to a nonsense value so we can later check in the table
  // element was filled in.
  for (int i=0; i<Path.NumSpecies(); i++)
    for (int j=0; j<Path.NumSpecies(); j++)
      PairMatrix(i,j) = (PairActionFitClass*)NULL;
  // Read pair actions files
  IOSectionClass PAIO;

  for (int i=0; i<numPairActions; i++) {
    // Allow for tilde-expansion in these files
    string name = ExpandFileName(PAFiles(i));
    assert(PAIO.OpenFile (name));
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
  if (HaveLongRange()) {
    assert (in.ReadVar("UseBackground", LongRange.UseBackground));
    LongRangePot.UseBackground = LongRange.UseBackground;
    LongRangeRPA.UseBackground = LongRange.UseBackground;
  }

//   if (longRange){
//     LongRange.Init(in);
//     if (UseRPA)
//       LongRangeRPA.Init(in);
//   }

  // Create nodal action objects
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

  ///Reading in information for David long range action
  if (PathData.Path.DavidLongRange){
    DavidLongRange.Read(in);
  }
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
      NodalActions(species) = 
	new FreeNodalActionClass (PathData, species);
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
    in.CloseSection();
  }
}



void
ActionsClass::Energy (double& kinetic, double &dUShort, double &dULong, 
		      double &node, double &vShort, double &vLong)
{
  int M = PathData.Path.NumTimeSlices()-1;
  kinetic = Kinetic.d_dBeta (0, M, 0);
  dUShort = ShortRange.d_dBeta (0, M, 0);
  dULong=0.0;
  if (PathData.Path.LongRange){
    if (UseRPA)
      dULong = LongRangeRPA.d_dBeta (0, M, 0);
    else
      dULong = LongRange.d_dBeta (0, M, 0);
  }
  node = 0.0;
  for (int species=0; species<PathData.Path.NumSpecies(); species++)
    if (NodalActions(species) != NULL)
      node += NodalActions(species)->d_dBeta(0, M, 0);

  vShort  = 0.0; vLong   = 0.0;  
  for (int slice=0; slice <= M; slice++) {
    double factor = ((slice==0)||(slice==M)) ? 0.5 : 1.0;
    vShort += factor * ShortRangePot.V(slice);
    vLong  += factor *  LongRangePot.V(slice);
  }
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
  ShortRangeApproximate.ShiftData(slicesToShift);
  LongRange.ShiftData(slicesToShift);
  LongRangeRPA.ShiftData(slicesToShift);
  DavidLongRange.ShiftData(slicesToShift);
  TIP5PWater.ShiftData(slicesToShift);
  ST2Water.ShiftData(slicesToShift);
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i)!=NULL)
      NodalActions(i)->ShiftData(slicesToShift);
}


void 
ActionsClass::AcceptCopy (int startSlice, int endSlice,
			       const Array<int,1> &activeParticles)
{
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->AcceptCopy (startSlice, endSlice);
}


void 
ActionsClass::RejectCopy (int startSlice, int endSlice,
			       const Array<int,1> &activeParticles)
{
  for (int i=0; i<NodalActions.size(); i++)
    if (NodalActions(i) != NULL)
      NodalActions(i)->RejectCopy (startSlice, endSlice);
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
  return longRange;
}
