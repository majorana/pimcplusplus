#include "ActionsClass.h"

#include "../PathDataClass.h"



#include "../Common/Ewald/OptimizedBreakup.h"
#include "../Common/Integration/GKIntegration.h"

///Actionsclass. Stores all the actsion
void ActionsClass::Read(IOSectionClass &in)
{ 
  PathClass &Path=PathData.Path;
  assert(in.ReadVar ("tau", Path.tau));
  assert(in.ReadVar ("MaxLevels", MaxLevels));
  cerr << "MaxLevels = " << MaxLevels << endl;

  if (!in.ReadVar ("UseRPA", UseRPA))
    UseRPA = false;
  if (UseRPA) 
    cerr << "Using RPA for long range action.\n";
  else      
    cerr << "Not using RPA for long range action.\n";


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
  bool longRange = false;
  for (int i=0; i<numPairActions; i++) {
    assert(PAIO.OpenFile (PAFiles(i)));
    PairArray(i) = ReadPAFit (PAIO, Path.tau, MaxLevels);
    longRange |= PairArray(i)->IsLongRange();
    bool paUsed=false;
    for (int spec1=0;spec1<Path.NumSpecies();spec1++)
      if (Path.Species(spec1).Type==PairArray(i)->Particle1.Name)
	for (int spec2=spec1;spec2<Path.NumSpecies();spec2++) 
	  if (Path.Species(spec2).Type==PairArray(i)->Particle2.Name) {
	    if (PairMatrix(spec1,spec2) != NULL) {
	      cerr << "More than one pair action for species types (" 
		   << PairArray(i)->Particle1.Name << ", "
		   << PairArray(i)->Particle2.Name << ")." << endl;
	      exit(-1);
	    }
	    cerr << "Found PAfile for pair (" 
		 << Path.Species(spec1).Name << ", "
		 << Path.Species(spec2).Name << ")\n";
	    PairMatrix(spec1,spec2) = PairArray(i);
	    PairMatrix(spec2,spec1) = PairArray(i);
	    paUsed = true;
	  }
    if (!paUsed) {
      cerr << "Warning:  Pair action for species types (" 
	   << PairArray(i)->Particle1.Name << ", "
 	   << PairArray(i)->Particle1.Name << ") not used.\n";
    }
    PAIO.CloseFile();
  }

  if (longRange){
    LongRange.Init(in);
    if (UseRPA)
      LongRangeRPA.Init(in);
  }
   
  // Now check to make sure all PairActions that we need are defined.
  for (int species1=0; species1<Path.NumSpecies(); species1++)
    for (int species2=0; species2<Path.NumSpecies(); species2++)
      if (PairMatrix(species1,species2) == NULL) {
	if ((species1 != species2) || 
	    (Path.Species(species1).NumParticles > 1)) {
	  cerr << "We're missing a PairAction for species1 = "
	       << Path.Species(species1).Name << " and species2 = "
	       << Path.Species(species2).Name << endl;
	  exit(1);
	}
      }

  // Create nodal action objects
  NodalActions.resize(PathData.Path.NumSpecies());
  for (int spIndex=0; spIndex<PathData.Path.NumSpecies(); spIndex++) {
    SpeciesClass &species = PathData.Path.Species(spIndex);
    if (species.GetParticleType() == FERMION) {
      if (species.NodeType == "FREE") 
	NodalActions (spIndex) = new FreeNodalActionClass(PathData, spIndex);
      else {
	cerr << "Unrecognized node type " << species.NodeType << ".\n";
	exit(EXIT_FAILURE);
      }
      NodalActions(spIndex)->Read(in);
    }
    else
      NodalActions(spIndex) = NULL;
  }

//   // Now create nodal actions for Fermions
//   NodalActions.resize(PathData.Path.NumSpecies());
//   for (int species=0; species<PathData.Path.NumSpecies(); species++) 
//     if (PathData.Path.Species(species).GetParticleType() == FERMION)
//       NodalActions(species) = new FPNodalActionClass(PathData, species);
//     else
//       NodalActions(species) = NULL;
  
  cerr << "Finished reading the action.\n"; 


  ///Reading in information for David long range action
  if (PathData.Path.DavidLongRange){
    DavidLongRange.Read(in);
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
