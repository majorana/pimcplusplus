#include "SpeciesClass.h"

bool SpeciesClass::Read(InputSectionClass &inSection)
{
  assert(inSection.ReadVar("Name",Name));
  assert(inSection.ReadVar("lambda",lambda));
  assert(inSection.ReadVar("NumParticles",NumParticles));
  assert(inSection.ReadVar("NumDim",NumDim));
  return true;
}
bool FermionClass::Read(InputSectionClass &inSection)
{
  bool success = SpeciesClass::Read(inSection);
  return success;
}


SpeciesClass* ReadSpecies(InputSectionClass &inSection)
{
  
  string typeString;
  SpeciesClass *mySpecies;
  assert(inSection.ReadVar("Type",typeString));
  if (typeString=="FERMION")
    mySpecies=new FermionClass();
  //  else if (typeString=="BOSON"){
  //    mySpecies=new BosonClass();
  //  }
  //  else if (typeString=="BOLTZMANON"){
  //    mySpecies=new BoltzmanonClass();
  //  }
  else {
    cerr<<"Species Type Unknown "<<typeString<<endl;
    exit(1);
  }
  assert(mySpecies->Read(inSection));
  return mySpecies;



}


double FermionClass::NodeAction(int ptcl,int LinkNum)
{


  return 0.0;
}




