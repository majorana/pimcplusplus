#include "PathClass.h"




void PathClass::Read (InputSectionClass &inSection)
{
  SetMode (BOTHMODE);
  double tau;
  assert(inSection.FindVar ("NumTimeSlices", TimeSliceNumber));
  assert(inSection.FindVar ("tau", tau));
  Array<double,1> tempBox;
  assert(inSection.FindVar ("Box",tempBox));
  assert(tempBox.size()==NDIM);
  for (int counter=0;counter<tempBox.size();counter++){
    Box(counter)=tempBox(counter);
  }
  assert(inSection.OpenSection("Particles"));
  int NumSpecies = inSection.CountSections ("Species");
  SpeciesArray.resize(NumSpecies);
  // First loop over species and read info about species
  for (int Species=0; Species < NumSpecies; Species++)
    {
      inSection.OpenSection("Species", Species);
      SpeciesClass *newSpecies = ReadSpecies (inSection);
      inSection.CloseSection(); // "Species"
      AddSpecies (newSpecies);
    }
  // Now actually allocate the path
  Allocate();
  // Now initilize the Path
  for (int speciesIndex=0; speciesIndex<NumSpecies; speciesIndex++);
  {
    SpeciesClass &species = *SpeciesArray(speciesIndex);
    assert(inSection.OpenSection("Species",Species));
    string InitPaths;
    inSection.ReadVar ("InitPaths", InitPaths);
    if (InitPaths == "RANDOM") {
      cerr << "Don't know how to do RANDOM yet.\n";
      exit(1);
    }
    else if (InitPaths == "FIXED") {
      Array<double,2> Positions;
      assert (inSection.ReadVar ("Positions", Positions));
      
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      for (int ptcl=species.FirtPtcl; 
	   ptcl<=species.LastPtcl; ptcl++)
	for (int slice=0; slice<NumTimeSlices; slice++) {
	  dVec pos;
	  pos = 0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
	  Path.SetPos(ptcl,slice,pos);
	}
    }
    else {
      cerr << "Unrecognize initialization strategy " 
	   << InitPaths << endl;
    }
  }

  CloseSection(); // "Particles"
  
}






