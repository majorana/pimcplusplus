#include "ActionsClass.h"

#include "../PathDataClass.h"



#include "../Common/Ewald/OptimizedBreakup.h"
#include "../Common/Integration/GKIntegration.h"

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

  PairMatrix = (PairActionFitClass*)NULL;
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
	for (int spec2=0;spec2<Path.NumSpecies();spec2++) 
	  if (Path.Species(spec2).Type==PairArray(i)->Particle2.Name) {
	    if (PairMatrix(spec1,spec2) != NULL) {
	      cerr << "More than one pair action for species types (" 
		   << PairArray(i)->Particle1.Name << ", "
		   << PairArray(i)->Particle1.Name << ")." << endl;
	      exit(-1);
	    }
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


