#include "PathClass.h"




void PathClass::Read (IOSectionClass &inSection)
{
  SetMode (BOTHMODE);
  double tau;
  assert(inSection.ReadVar ("NumTimeSlices", TotalNumSlices));
  assert(inSection.ReadVar ("tau", tau));
  Array<double,1> tempBox;
  Array<bool,1> tempPeriodic;
  TinyVector<bool,NDIM> periodic;
  assert (inSection.ReadVar ("IsPeriodic", tempPeriodic));
  assert (tempPeriodic.size() == NDIM);
  bool needBox = false;
  for (int i=0; i<NDIM; i++) {
    needBox = needBox || tempPeriodic(i);
    periodic(i) = tempPeriodic(i);
  }
  SetPeriodic (periodic);
  if (needBox) {
    assert(inSection.ReadVar ("Box", tempBox));
    cerr << "Using periodic boundary conditions.\n";
    assert(tempBox.size()==NDIM);
    for (int counter=0;counter<tempBox.size();counter++)
      Box(counter)=tempBox(counter);
    SetBox (Box);
    kCut=-1;
    inSection.ReadVar("kCutoff",kCut);
  }
  else 
    cerr << "Using free boundary conditions.\n";

  assert(inSection.OpenSection("Particles"));
  int NumSpecies = inSection.CountSections ("Species");
  cerr<<"we have this many sections: "<<NumSpecies<<endl;
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
  // Now initialize the Path
  for (int speciesIndex=0; speciesIndex<NumSpecies; speciesIndex++)
  {
    SpeciesClass &species = *SpeciesArray(speciesIndex);
    assert(inSection.OpenSection("Species", speciesIndex));
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
      for (int ptcl=species.FirstPtcl; 
	   ptcl<=species.LastPtcl; ptcl++)
	for (int slice=0; slice<NumTimeSlices(); slice++) {
	  dVec pos;
	  pos = 0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
	  SetPos(slice,ptcl,pos);
	}      
    }
    else {
      cerr << "Unrecognize initialization strategy " 
	   << InitPaths << endl;
    }
    inSection.CloseSection();
  }

  inSection.CloseSection(); // "Particles"
  
}

PathClass::SetupkVectors()
{
 
  int currVec=0;
  for (kx=-ceil(kCut/kBox(0));kx<ceil(kCut/kBox(0));kx++){
    for (ky=-ceil(kCut/kBox(1));ky<ceil(kCut/kBox(1));ky++){
        for (kz=-ceil(kCut/kBox(1));kx<ceil(kCut/kBox(2));kz++){
	  kMag2=kx*kx+ky*ky+kz*kz;
	  if (kMag2<kCut*kCut && (kx>0 || (kx==0 & ky>0) || (kx==0 && ky==0 && kz>0))){
	    kVectors(currVec)(0)=kx;
	    kVectors(currVec)(1)=ky;
	    kVectors(currVec)(2)=kz;
	  }
	}
    }
  }      
}




