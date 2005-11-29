#include "Forces.h"

void
ForcesClass::SetSpecies(int speciesNum)
{
  TimesCalled = 0;
  SpeciesNum = speciesNum;
  SpeciesClass &species = PathData.Path.Species(speciesNum);
  Forces.resize(species.NumParticles);
  dVec zero;
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  Forces = zero;
  Ptcls.resize(species.NumParticles);
  for (int i=0; i<Ptcls.size(); i++)
    Ptcls(i) = i + species.FirstPtcl;
  ForcesArray.resize(Ptcls.size(),NDIM);
}


void
ForcesClass::Accumulate()
{
  TimesCalled++;

  if ((TimesCalled % DumpFreq) == 0)
    WriteBlock();

  if ((TimesCalled % Freq) == 0) {
//     Array<dVec,1> Fanalytic(Ptcls.size()), FFD(Ptcls.size());
//     dVec zero;
//     for (int i=0; i<NDIM; i++) zero[i] = 0.0;
//     Fanalytic = zero;
//     FFD = zero;
//     PathData.Actions.GetForces(Ptcls, Fanalytic);
//     PathData.Actions.GetForcesFD(Ptcls, FFD);
//     cerr << "Forces:\n";
//     for (int i=0; i<Ptcls.size(); i++) {
//       for (int dim=0; dim<NDIM; dim++)
// 	fprintf (stderr, "%12.6f ", Fanalytic(i)[dim]);
//       for (int dim=0; dim<NDIM; dim++)
// 	fprintf (stderr, "%12.6f ", FFD(i)[dim]);
//       fprintf (stderr, "\n");
//     }
    PathData.Actions.GetForces(Ptcls, Forces);
    Counts++;
  }
}


void
ForcesClass::WriteBlock()
{
  double norm = 1.0/(double)Counts;
  for (int i=0; i<Forces.size(); i++)
    for (int j=0; j<NDIM; j++) 
      ForcesArray(i,j) = norm * Forces(i)[j];
  ForcesVar.Write(ForcesArray);
  dVec zero;
  for (int i=0; i<NDIM; i++)
    zero[i] = 0.0;
  Forces = zero;
  Counts = 0;
}


void
ForcesClass::WriteInfo()
{
  IOSection.WriteVar("SpeciesNum", SpeciesNum);
  Array<string,1> actionTypes(2);
  actionTypes(0) = "ShortRange";
  actionTypes(1) = "LongRange";
  IOSection.WriteVar("ActionTypes", actionTypes);
}


void
ForcesClass::Read(IOSectionClass &in)
{
  string speciesString;
  assert (in.ReadVar("Species", speciesString));
  SpeciesNum = PathData.Path.SpeciesNum(speciesString);
  SetSpecies(SpeciesNum);
  assert (in.ReadVar("freq", Freq));
  assert (in.ReadVar("dumpFreq", DumpFreq));
  assert (in.ReadVar("Name", Name));
  WriteInfo();
}
