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

  in.ReadVar ("UseRPA", UseRPA);
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
    // cerr << "i = " << i << endl;
    assert(PAIO.OpenFile (PAFiles(i)));
    PairArray(i) = ReadPAFit (PAIO, Path.tau, MaxLevels);
    longRange |= PairArray(i)->IsLongRange();
    int type1 = Path.SpeciesNum(PairArray(i)->Particle1.Name);
    int type2 = Path.SpeciesNum(PairArray(i)->Particle2.Name);
    if (type1==-1) {
      cerr << "Unrecognized type \""
	   << PairArray(i)->Particle1.Name << "\".\n";
      abort();
    }
    if (type2==-1) {
      cerr << "Unrecognized type \""
	   << PairArray(i)->Particle2.Name << "\".\n";
      abort();
    }
    PairMatrix(type1,type2) = PairArray(i);
    PairMatrix(type2,type1) = PairArray(i);
    PAIO.CloseFile();
  }

  if (longRange){
    LongRange.Init(in);
    if (UseRPA)
      LongRange.Init(in);
  }
 

//     if (UseRPA) {
//       cerr << "Doing RPA correction...\n";
//       SetupRPA();
//       cerr << "done.\n";
//     }
//     // Print out some debug info
//     for (int i=0; i<numPairActions; i++) {
//       string fname = PairArray(i)->Particle1.Name + "-" +
// 	PairArray(i)->Particle2.Name + ".dat";
//       FILE *fout = fopen (fname.c_str(), "w");
//       for (double q=0.0; q<20000.0; q+=1.0) {
// 	double Udiag = PairArray(i)->Udiag(q, 0);
// 	double Ulong = PairArray(i)->Ulong(0)(q);
// 	double dUdiag = PairArray(i)->dUdiag(q, 0);
// 	double dUlong = PairArray(i)->dUlong(0)(q);
// 	double V = PairArray(i)->V(q);
// 	double Vlong = PairArray(i)->Vlong(q);
// 	fprintf (fout, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", 
// 		 q, Udiag, Ulong, dUdiag, dUlong, V, Vlong);
//       }
//       fclose (fout);
//     }
//  }
  


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

//   // Now create nodal actions for Fermions
//   NodalActions.resize(PathData.Path.NumSpecies());
//   for (int species=0; species<PathData.Path.NumSpecies(); species++) 
//     if (PathData.Path.Species(species).GetParticleType() == FERMION)
//       NodalActions(species) = new FPNodalActionClass(PathData, species);
//     else
//       NodalActions(species) = NULL;

  cerr << "Finished reading the action.\n"; 
}


