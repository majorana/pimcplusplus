#include "PathClass.h"

void SetupPath (PathClass &path)
{
  SetMode (NEWMODE);
  path.TotalNumSlices = 10;
  TinyVector<bool,3> periodic;
  periodic = true, true, true;
  path.SetPeriodic (periodic);
  dVec box;
  box = 1.0, 1.2, 1.4;
  path.SetBox(box);
  
  path.kCutoff = 30.0;
  
  FermionClass *protons=new FermionClass();
  FermionClass *electrons=new FermionClass();
  protons->Name = "protons";
  protons->lambda = 0.0;
  protons->NumParticles = 4;
  protons->NumDim = 3;

  electrons->Name = "electrons";
  electrons->lambda = 0.5;
  electrons->NumParticles = 4;
  electrons->NumDim = 3;

  path.AddSpecies (protons);
  path.AddSpecies (electrons);

  path.Allocate();

  /// Initialize the path positions for NaCl structure.
  for (int slice=0; slice <= path.TotalNumSlices; slice++) {
    path(slice, 0) = 0.0, 0.5, 0.5;
    path(slice, 1) = 0.5, 0.0, 0.5;
    path(slice, 2) = 0.0, 0.0, 0.0;
    path(slice, 3) = 0.5, 0.5, 0.0;
    path(slice, 4) = 0.0, 0.5, 0.0;
    path(slice, 5) = 0.5, 0.0, 0.0;
    path(slice, 6) = 0.0, 0.0, 0.5;
    path(slice, 7) = 0.5, 0.5, 0.5;
  }

  path.Path.AcceptCopy();
}


void TestRho_k(PathClass &path)
{
  SetMode(NEWMODE);
  for (int slice=0; slice<path.NumTimeSlices(); slice++) 
    for (int species=0; species<path.NumSpecies(); species++)
      path.CalcRho_ks_Fast(slice,species);
  SetMode(OLDMODE);
  for (int slice=0; slice<path.NumTimeSlices(); slice++) 
    for (int species=0; species<path.NumSpecies(); species++)
      path.CalcRho_ks_Slow(slice,species);

  for (int slice=0; slice<path.NumTimeSlices(); slice++) 
    for (int species=0; species<path.NumSpecies(); species++)
      for (int ki=0; ki<path.kVecs.size(); ki++)
	cerr << "Fast = " << path.Rho_k[0](slice, species, ki)
	     << " Slow = " << path.Rho_k[1](slice, species, ki) << endl;
}


main()
{
  PIMCCommunicatorClass comm;
  PathClass path(comm);
  SetupPath(path);
  TestRho_k(path);
}
