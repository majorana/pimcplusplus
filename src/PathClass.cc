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
  }
  else 
    cerr << "Using free boundary conditions.\n";

  // Read in the k-space radius.  If we don't have that,
  // we're not long-ranged.
  LongRange = inSection.ReadVar("kCutoff", kCutoff);

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


inline bool Include(dVec k)
{
  assert (NDIM == 3);
  if (k[0] > 0.0)
    return true;
  else if ((k[0]==0.0) && (k[1]>0.0))
    return true;
  else if ((k[0]==0.0) && (k[1]==0.0) && (k[2] > 0.0))
    return true;
  else
    return false;
}
 

void PathClass::Allocate()
{
  assert(TotalNumSlices>0);
  int myProc=Communicator.MyProc();
  int numProcs=Communicator.NumProcs();
  ///Everybody gets the same number of time slices if possible.
  ///Otherwise the earlier processors get the extra one slice 
  ///until we run out of extra slices.
  ///The last slice on processor i is the first slices on processor i+1
  MyNumSlices=TotalNumSlices/numProcs+1+(myProc<(TotalNumSlices % numProcs));
  cerr<<"Numprocs is "<<numProcs<<endl;
  cerr<<"mynumslices: "<<MyNumSlices<<endl;
  
  
  int numParticles = 0;
  /// Set the particle range for the new species
  for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
    SpeciesArray(speciesNum)->FirstPtcl = numParticles;
    numParticles=numParticles + SpeciesArray(speciesNum)->NumParticles;
    SpeciesArray(speciesNum)->LastPtcl= numParticles-1;
  }
  Path.Resize(MyNumSlices,numParticles);
  Permutation.Resize(numParticles);
  SpeciesNumber.resize(numParticles);
  DoPtcl.resize(numParticles);
  /// Assign the species number to the SpeciesNumber array
  for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
    for (int i=SpeciesArray(speciesNum)->FirstPtcl; 
	 i<= SpeciesArray(speciesNum)->LastPtcl; i++)
      SpeciesNumber(i) = speciesNum;
  }
  //Sets to the identity permutaiton 
  for (int ptcl=0;ptcl<Permutation.NumParticles();ptcl++){
    Permutation.Set(ptcl,ptcl);
  }

  SetupkVecs();
  Rho_k.resize(NumSpecies(), MyNumSlices, kVecs.size());
}


void PathClass::SetupkVecs()
{
  dVec kBox;
  for (int i=0; i<NDIM; i++)
    kBox[i] = 2.0*M_PI/Box(i);
  
  int numVecs=0;

  assert (NDIM == 3);
  int nx = (int) ceil(1.25*kCutoff / kBox[0]);
  int ny = (int) ceil(1.25*kCutoff / kBox[1]);
  int nz = (int) ceil(1.25*kCutoff / kBox[2]);
  
  dVec k;
  for (int ix=-nx; ix<=nx; ix++) {
    k[0] = ix*kBox[0];
    for (int iy=-nx; iy<=ny; iy++) {
      k[1] = iy*kBox[1];
      for (int iz=-nx; iz<=nz; iz++) {
	k[2] = iz*kBox[2];
	if ((dot(k,k)<kCutoff*kCutoff) && Include(k))
	  numVecs++;
      }
    }
  }
  kVecs.resize(numVecs);
  numVecs = 0;
  for (int ix=-nx; ix<=nx; ix++) {
    k[0] = ix*kBox[0];
    for (int iy=-nx; iy<=ny; iy++) {
      k[1] = iy*kBox[1];
      for (int iz=-nx; iz<=nz; iz++) {
	k[2] = iz*kBox[2];
	if ((dot(k,k)<kCutoff*kCutoff) && Include(k)) {
	  kVecs(numVecs) = k;
	  numVecs++;
	}
      }
    }
  }
}


void PathClass::CalcRho_ks()
{
  for (int speciesIndex=0; speciesIndex<NumSpecies(); speciesIndex++)
    for (int slice=0; slice<MyNumSlices; slice++)
      for (int ki=0; ki<kVecs.size(); ki++) {
	complex<double> rho;
	rho = 0.0;
	for (int ptcl=Species(speciesIndex).FirstPtcl; 
	     ptcl <= Species(speciesIndex).LastPtcl; ptcl++) {
	  const dVec &r = (*this)(slice, ptcl);
	  double phase = dot(r, kVecs(ki));
	  rho += complex<double> (cos(phase), sin(phase));
	}
	Rho_k(speciesIndex, slice, ki) = rho;
      }
}
	  
	


