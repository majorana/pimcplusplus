#include "Common/IO/InputOutput.h"



int main()
{
  OutputSectionHDF5Class outSection;
  outSection.OpenFile("hydrogen.h5");
  outSection.WriteVar("tau",1.0);
  outSection.OpenSection("System");
  outSection.WriteVar("NumTimeSlices",50);
  Array <double,1> Box(3);
  Box=5.0,5.0,5.0;
  outSection.WriteVar("Box",Box);
  outSection.OpenSection("Particles");

  outSection.OpenSection("Species");
  outSection.WriteVar("Name","electron");
  outSection.WriteVar("lambda",0.5);
  outSection.WriteVar("Type","FERMION");
  outSection.WriteVar("NumParticles",1);
  outSection.WriteVar("NumDim",3);
  outSection.WriteVar("InitPaths","FIXED");
  Array <double,2> Positions(1,3);
  Positions= 0.1, 0.0, 0.0;
  outSection.WriteVar("Positions",Positions);
  outSection.CloseSection();

  outSection.OpenSection("Species");
  outSection.WriteVar("Name","proton");
  outSection.WriteVar("lambda",0.0);
  outSection.WriteVar("Type","FERMION");
  outSection.WriteVar("NumParticles",1);
  outSection.WriteVar("NumDim",3);
  outSection.WriteVar("InitPaths","FIXED");
  Positions= 0.0, 0.0, 0.0;
  outSection.WriteVar("Positions",Positions);
  outSection.CloseSection();

  outSection.CloseSection(); //Particles;
  outSection.CloseSection(); // System;
  outSection.OpenSection("Moves");

  outSection.OpenSection("Move");
  outSection.WriteVar("type","Bisection");
  Array<string,1> ActiveSpecies(1);
  ActiveSpecies="electron";
  outSection.WriteVar("ActiveSpecies",ActiveSpecies);
  outSection.WriteVar("NumParticlesToMove",1);
  outSection.WriteVar("NumLevels",3);
  outSection.CloseSection(); //Move

  outSection.CloseSection(); //Moves

  outSection.OpenSection("Action");
  outSection.OpenSection("PairAction");
  outSection.WriteVar("type1","electron");
  outSection.WriteVar("type2","proton");
  outSection.WriteVar("dmfile","../inputs/ep_beta1.0.dm");
  outSection.CloseSection();//PairAction

  outSection.CloseSection(); //Action;
  
  outSection.CloseFile();


   

}
