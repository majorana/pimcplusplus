#include "ActionClass.h"
#include "PathDataClass.h"

void SetupPath (PathClass &path)
{
  SetMode (NEWMODE);
  path.TotalNumSlices = 10;
  TinyVector<bool,3> periodic;
  periodic = true, true, true;
  path.SetPeriodic (periodic);
  dVec box;
  box = 2.0, 2.0, 2.0;
  path.SetBox(box);
  
  path.kCutoff = 25.0;
  
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
    path(slice, 0) = 0.0, 1.0, 1.0;
    path(slice, 1) = 1.0, 0.0, 1.0;
    path(slice, 2) = 0.0, 0.0, 0.0;
    path(slice, 3) = 1.0, 1.0, 0.0;
    path(slice, 4) = 0.0, 1.0, 0.0;
    path(slice, 5) = 1.0, 0.0, 0.0;
    path(slice, 6) = 0.0, 0.0, 1.0;
    path(slice, 7) = 1.0, 1.0, 1.0;
  }

  path.Path.AcceptCopy();
}


void TestRho_k(PathClass &path)
{
  #include <time.h>
  clock_t start, end;
  double deltat;

  SetMode(NEWMODE);
  start = clock();
  for (int i=0; i<10; i++) 
    for (int slice=0; slice<path.NumTimeSlices(); slice++) 
      for (int species=0; species<path.NumSpecies(); species++)
	path.CalcRho_ks_Fast(slice,species);
  end = clock();
  deltat = (double)(end-start)/1.0e6;
  cerr << "Fast speed = " << 1.0e1/deltat << endl;
  SetMode(OLDMODE);
  start = clock();
  for (int i=0; i<10; i++) 
    for (int slice=0; slice<path.NumTimeSlices(); slice++) 
      for (int species=0; species<path.NumSpecies(); species++)
	path.CalcRho_ks_Slow(slice,species);
  end = clock();
  deltat = (double)(end-start)/1.0e6;
  cerr << "Slow speed = " << 1.0e1/deltat << endl;
  for (int slice=0; slice<path.NumTimeSlices(); slice++) 
    for (int species=0; species<path.NumSpecies(); species++)
      for (int ki=0; ki<path.kVecs.size(); ki++)
	if (fabs(path.Rho_k[0](slice,species,ki).real() - 
		 path.Rho_k[1](slice,species,ki).real())> 1.0e-13 ||
	    fabs(path.Rho_k[0](slice,species,ki).imag() - 
		 path.Rho_k[1](slice,species,ki).imag())> 1.0e-13)
	  cerr << "Fast/slow discrepancy in TestRho_k\n";
  // Make copies exactly the same.
  path.Rho_k[0] = path.Rho_k[1];
}

void SetupAction(ActionClass &action,PathDataClass &pathData)
{
  IOSectionClass in;
  in.OpenFile ("EwaldTest.in");
  action.tau = 1.0;
  action.MaxLevels=3;
  Array<string,1> PAFiles(3);
  PAFiles="p-p.PairAction","e-e.PairAction","p-e.PairAction";
  int numPairActions=3;
  action.PairActionVector.resize(3);
  IOSectionClass PAIO;
  for (int i=0;i<numPairActions;i++){
    action.PairActionVector(i) = new PAclassicalFitClass;
    assert(PAIO.OpenFile(PAFiles(i)));
    action.PairActionVector(i)=ReadPAFit(PAIO,action.tau,action.MaxLevels);
    cerr<<"Kvecs size: "<<pathData.Path.kVecs.size()<<endl;
    action.PairActionVector(i)->DoBreakup(pathData.Path.GetBox(),
					  pathData.Path.kVecs);
    PAIO.CloseFile();
  }
  action.PairMatrix.resize(2,2);  
  action.PairMatrix(0,0)=0;
  action.PairMatrix(1,1)=1;
  action.PairMatrix(0,1)=2;
  action.PairMatrix(1,0)=2;

  
}


void MadelungTest(ActionClass &action)
{
  double longRange=action.CalcLRAction(0,0);
  Array<int,1> changedParticles(8);
  changedParticles=0,1,2,3,4,5,6,7;
  double shortRange=action.calcTotalAction(0,1,changedParticles,0);
  cerr<<"LongRange, ShortRange, Total: "<<longRange<<", "<<shortRange;
  cerr<<", "<<longRange+shortRange<<endl;
}

main()
{
  PathDataClass pathData;
  SetupPath(pathData.Path);
  ActionClass action(pathData);
  SetupAction (action,pathData);
  TestRho_k(pathData.Path);
  MadelungTest(action);

}
