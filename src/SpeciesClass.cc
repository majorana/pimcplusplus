#include "SpeciesClass.h"

bool SpeciesClass::Read(IOSectionClass &inSection)
{
  assert(inSection.ReadVar("Name",Name));
  assert(inSection.ReadVar("lambda",lambda));
  assert(inSection.ReadVar("NumParticles",NumParticles));
  assert(inSection.ReadVar("NumDim",NumDim));
  assert(inSection.ReadVar("Type",Type));
  return true;
}

bool FermionClass::Read(IOSectionClass &inSection)
{
  bool success = SpeciesClass::Read(inSection);
  assert (inSection.ReadVar ("NodeType", NodeType));
  return success;
}


SpeciesClass* ReadSpecies(IOSectionClass &inSection)
{
  
  string statisticsString;
  SpeciesClass *mySpecies;
  assert(inSection.ReadVar("Statistics",statisticsString));
  //  cerr<<"My type string is "<<typeString<<endl;
  if (statisticsString=="FERMION")
    mySpecies=new FermionClass();
  else if (statisticsString=="BOSON"){
    mySpecies=new BosonClass();
  }
  else if (statisticsString=="BOLTZMANNON"){
    mySpecies=new BoltzmannonClass();
  }
  else {
    cerr<<"Species Statistics Unknown "<<statisticsString<<endl;
    exit(1);
  }
  assert(mySpecies->Read(inSection));
  return mySpecies;



}





