#include "SpeciesClass.h"

bool SpeciesClass::Read(IOSectionClass &inSection)
{
  assert(inSection.ReadVar("Name",Name));
  assert(inSection.ReadVar("lambda",lambda));
  assert(inSection.ReadVar("NumParticles",NumParticles));
  assert(inSection.ReadVar("NumDim",NumDim));
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
  
  string typeString;
  SpeciesClass *mySpecies;
  assert(inSection.ReadVar("Type",typeString));
  //  cerr<<"My type string is "<<typeString<<endl;
  if (typeString=="FERMION")
    mySpecies=new FermionClass();
  else if (typeString=="BOSON"){
    mySpecies=new BosonClass();
  }
  else if (typeString=="BOLTZMANNON"){
    mySpecies=new BoltzmannonClass();
  }
  else {
    cerr<<"Species Type Unknown "<<typeString<<endl;
    exit(1);
  }
  assert(mySpecies->Read(inSection));
  return mySpecies;



}





