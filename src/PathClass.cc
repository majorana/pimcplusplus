#include "PathClass.h"


void PathClass::Read(InputSectionClass *section)
{
  int tempTimeSlices;
  int specNum=-1;
  Array<string,1> initPaths;
  initPaths.resize(10); ////HACK! CAN't BE MORE THEN 10 PARTICLE TYPES
  if (!(section->ReadVar("NumTimeSlices",tempTimeSlices))){
    cerr<<"Error reading the Time Slice Number!!!\n";
  }
  SetTimeSlices(tempTimeSlices);
  if (!(section->FindSection("Species",section,false))){
    cerr<<"Error finding the species section!!!!\n";
  }
  section->Rewind();
  cerr<<"Just rewound\n";
  InputSectionClass* particleSection;
  while (section->FindSection("Particle",particleSection,false)){
    cerr<<"I've found a section!!!"<<endl;
    specNum++;
    string typeString;
    particleSection->ReadVar("type",typeString);
    cerr<<"my type string is "<<typeString<<endl;
    if (typeString=="fermion"){
      FermionClass *myFermionPtr = new FermionClass();
      if (!(particleSection->ReadVar("lambda",(*myFermionPtr).lambda))){
	cerr<<"Error reading the lambda\n";
      }
      if (!(particleSection->ReadVar("NumParticles",(*myFermionPtr).NumParticles))){
	cerr<<"Error reading the number of particles\n";
      }
      cerr<<"Particles are "<<(*myFermionPtr).NumParticles;
      if (!(particleSection->ReadVar("InitPath",initPaths(specNum)))){
	cerr<<"Error reading the path initialization\n";
      }     
      AddSpecies(myFermionPtr);
    }
    else {
      cerr<<"The type isn't fermion and I don't know how to deal with it\n";
    }
      cerr<<"and here 2\n";
  }
  cerr<<"and here 3\n";
  Allocate();
  SetMode(BOTHMODE);
  double tau;
  if (!(section->ReadVar("tau",tau))){
    cerr<<"Can't read tau\n";
  }
  double sigma=sqrt(2*0.5*tau);
  for (int counter3=0;counter3<SpeciesArray.size();counter3++){
    for (int counter=SpeciesArray(counter3)->FirstPtcl;
	 counter<=SpeciesArray(counter3)->LastPtcl;counter++){
      for (int counter2=0;counter2<NumTimeSlices();counter2++){
	if (initPaths(counter3)=="fixed"){ //I should probably swap loops here
	  dVec zeroVector=0;
	  SetPos(counter2,counter,zeroVector);
	}
	else if (initPaths(counter3)=="random"){
	  dVec fermionVector=GaussianRandomVec(sigma);	  
	  SetPos(counter2,counter,fermionVector);
	}
	else {
	  cerr<<"Warning. You haven't told me how to initialize the paths.";
	}
      }
    }
  }
    


}


