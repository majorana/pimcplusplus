#include "PathClass.h"
#include "Actions/ActionsClass.h"

void PathClass::ReadOld(string fileName,bool replicate)
{
// This was modified on Jan 19 2005 to read in a set of classical (P=1) configs and duplicate them to produce a set of PIMC (P>1) configs.  -jg
  IOSectionClass inFile;
  assert (inFile.OpenFile(fileName.c_str()));
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<double,4> oldPaths; //(58,2560,2,3);
  
  assert(inFile.ReadVar("Path",oldPaths));
  cerr << "My paths are of size"  << oldPaths.extent(0) << " "
       << oldPaths.extent(1)<<" " << oldPaths.extent(2) << endl;
  
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
	if (replicate){
          pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,0,dim);
        }
        else{
	  pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,slice,dim);
        }
      Path(slice,ptcl) = pos;
    }      
  }
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();
}

void PathClass::RefDistDisp (int slice, int refPtcl, int ptcl,
			     double &dist, dVec &disp)
{
  disp = Path(slice, ptcl)- RefPath(refPtcl);
  
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
    if (!(-Box(i)/2.0<=disp(i))){
      cerr<<"ERROR: "<<Box(i)<<" "<<disp(i)<<" "
	  <<slice<<" "<<ptcl<<" "<<refPtcl<<" "
	  <<BoxInv(i)<<Path(slice,ptcl)<<" "<<endl;
//      sleep(5000);
    }
    assert(-Box(i)/2.0<=disp(i));
    assert(disp(i)<=Box(i)/2.0);
  }
  dist = sqrt(dot(disp,disp));

#ifdef DEBUG
  dVec DBdisp = Path(slice, ptcl) -RefPath(refPtcl);
  for (int i=0; i<NDIM; i++) {
    while (DBdisp(i) > 0.5*Box(i))
      DBdisp(i) -= Box(i);
    while (DBdisp(i) < -0.5*Box(i)) 
      DBdisp(i) += Box(i);
    if (fabs(DBdisp(i)-disp(i)) > 1.0e-12){ 
      cerr<<DBdisp(i)<<" "<<disp(i)<<endl;
    }
    //    assert (fabs(DBdisp(i)-disp(i)) < 1.0e-12);
  }
#endif
}



/// Constructs a Levi flight beginning in the vec(0) and ending
/// in vec(N-1) if vec has length N.  This is a path which samples
/// exactly the free particle density matrix.
void PathClass::LeviFlight (Array<dVec,1> &vec, double lambda)
{
  // HACK HACK HACK
  //  tau *= 0.000001;

  int N = vec.size();
  for (int slice=1; slice<(N-1); slice++) {
    double di = (double)slice;
    double delta = (double)(N-slice-1);
    dVec center = (1.0/(delta+1.0))*(delta*vec(slice-1) + vec(N-1));
    double taueff = tau*(1.0 - 1.0/delta);
    double sigma = sqrt (2.0*lambda*taueff);
    Random.CommonGaussianVec(sigma, vec(slice));
    vec(slice) += center;
  }
  // DEBUG DEBUG DEBUG DEBUG DEBUG!!!!
  FILE *fout = fopen ("Levi.dat", "w");
  for (int i=0; i<N; i++) {
    for (int dim=0; dim<NDIM; dim++) 
      fprintf (fout, "%1.12e ", vec(i)[dim]);
    fprintf (fout, "\n");
  }
  fclose (fout);
}

void 
PathClass::NodeAvoidingLeviFlight (int speciesNum, Array<dVec,1> &R0)
{
  SpeciesClass &species = Species(speciesNum);
  double lambda = species.lambda;
  bool haveNodeAction = Actions.NodalActions(speciesNum)!=NULL;

  int myFirstSlice, myLastSlice, myProc;
  myProc = Communicator.MyProc();
  SliceRange (myProc, myFirstSlice, myLastSlice);

  int numPtcls = species.NumParticles;
  Array<dVec,1> prevSlice (numPtcls), newSlice(numPtcls);
  for (int ptcl=0; ptcl<numPtcls; ptcl++) {
    prevSlice (ptcl) = R0(ptcl);
    if (myProc == 0)
      (*this)(0, ptcl+species.FirstPtcl) = R0(ptcl);
    RefPath(ptcl+species.FirstPtcl) = R0(ptcl);
  }

  /// Check slice 0 -- important for ground state nodes
  bool swapFirst = false;
  if (myProc == 0)
    if (haveNodeAction)
      if (! Actions.NodalActions(speciesNum)->IsPositive(0)) {
	cerr << "Initially negative node action.  Swapping two particles "
	     << "for species " << species.Name << ".\n";
	swapFirst = true;
      }
  Communicator.Broadcast (0,swapFirst);

  if (swapFirst) {
    dVec tmp = R0(0);
    R0(0) = R0(1);
    R0(1) = tmp;
    RefPath(species.FirstPtcl)     = R0(0);
    RefPath(species.FirstPtcl+1)   = R0(1);
    prevSlice(0)                   = R0(0);
    prevSlice(1)                   = R0(1);
    if (myProc==0) {
      (*this)(0,species.FirstPtcl)   = R0(0);
      (*this)(0,species.FirstPtcl+1) = R0(1);  
    }
  }

  int N = TotalNumSlices+1;
  for (int slice=1; slice<N; slice++) {
    int sliceOwner = SliceOwner(slice);
    int myProc = Communicator.MyProc();
    int relSlice = slice-myFirstSlice;
    
    double delta = (double)(N-slice-1);
    double taueff = tau*(1.0 - 1.0/(delta+1.0));
    double sigma = sqrt (2.0*lambda*taueff);
    bool positive = false;
    
    do {
      // Randomly construct new slice
      for (int ptcl=0; ptcl<numPtcls; ptcl++) {
	dVec center = (1.0/(delta+1.0))*(delta*prevSlice(ptcl) + R0(ptcl));
	Random.CommonGaussianVec(sigma, newSlice(ptcl));
	newSlice(ptcl) += center;
      }
      
      // Now check the nodal sign if we're a fermion species
      if ((Actions.NodalActions(speciesNum) == NULL) || (slice == (N-1)))
	positive = true;
      else {
	if (sliceOwner == myProc ) {
	  // Now assign to Path
	  for (int ptcl=0; ptcl<numPtcls; ptcl++)
	    (*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
	  positive = 
	    Actions.NodalActions(speciesNum)->IsPositive(relSlice);
	}
	// Now broadcast whether or not I'm positive to everyone
	Communicator.Broadcast(sliceOwner, positive);
      }
  } while (!positive);
//       if (!positive){
// 	cerr << "Negative sign at slice " << slice << "\n";
// 	// Find two closest particles in species.
// 	double  minDist = 1.0e300;
// 	int min1, min2;
// 	for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) 
// 	  for (int ptcl2=ptcl1+1; ptcl2<numPtcls; ptcl2++) {
// 	    double dist;
// 	    dVec disp;
// 	    DistDisp(slice, ptcl1+species.FirstPtcl, ptcl2+species.LastPtcl, dist, disp);
// 	    if (dist < minDist) {
// 	      minDist = dist;
// 	      min1 = ptcl1; 
// 	      min2 = ptcl2;
// 	    }
// 	  }
// 	// Now, swap them
// 	cerr << "  Swapping ptcl #" << min1 << " and #" << min2 << endl;
// 	dVec tmp = newSlice(min1);
// 	newSlice(min1) = newSlice(min2);
// 	newSlice(min2) = tmp;
// 	// Swap in path if I'm the slice owner
//       }
//     }
    // Copy slice into Path if I'm the slice owner.
    if ((slice>=myFirstSlice) && (slice<=myLastSlice)) 
      for (int ptcl=0; ptcl<numPtcls; ptcl++) {
	(*this)(relSlice, ptcl+species.FirstPtcl) = newSlice(ptcl);
	//cerr << "setting slice " << relSlice 
	//     << " for species " << species.Name << endl;
      }

    // Check to make sure we're positive now.
    if ((Actions.NodalActions(speciesNum) != NULL) && (slice!=(N-1))) {
      bool positive;
      if (sliceOwner == myProc) {
	positive = Actions.NodalActions(speciesNum)->IsPositive(relSlice);
	if (!positive) {
	  cerr << "Still not postive at slice " << slice << ".\n";
	  abort();
	}
      }
    }

    // continue on to next slice
    prevSlice = newSlice;
  }
  Array<int,1> changedParticles(1);
  if (haveNodeAction) {
    double action = 
      Actions.NodalActions(speciesNum)->Action(0, N-1, 
					       changedParticles,0);
    cerr << "Nodal Action after Levi flight = " << action << endl;
  }

  char fname[100];
  snprintf (fname, 100, "%s.dat", species.Name.c_str());
  FILE *fout = fopen (fname, "w");
  for (int slice=0; slice<N; slice++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) 
      for (int i=0; i<NDIM; i++)
	fprintf (fout, "%1.12e ", (*this)(slice, ptcl)[i]);
    fprintf (fout, "\n");
  }
  fclose(fout);
  cerr << "My first particle = " << species.FirstPtcl << endl;
}

void PathClass::Read (IOSectionClass &inSection)
{
  SetMode (NEWMODE);
  Weight=1;
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
  if (!inSection.ReadVar("OpenLoops",OpenPaths))
    OpenPaths=false;

  // Read in the k-space radius.  If we don't have that,
  // we're not long-ranged.
  LongRange = inSection.ReadVar("kCutoff", kCutoff);
  if (!(inSection.ReadVar("DavidLongRange",DavidLongRange))){
    DavidLongRange=false;
  }
  if (DavidLongRange)
    cerr<<"I am doing DAVID LONG RANGE!"<<endl;
  assert(inSection.OpenSection("Particles"));
  int numSpecies = inSection.CountSections ("Species");
  cerr<<"we have this many sections: "<<numSpecies<<endl;
   // First loop over species and read info about species
  for (int Species=0; Species < numSpecies; Species++) {
    inSection.OpenSection("Species", Species);
    SpeciesClass *newSpecies = ReadSpecies (inSection);
    inSection.CloseSection(); // "Species"
    AddSpecies (newSpecies);
  }
  inSection.CloseSection(); // Particles
  // Now actually allocate the path
  Allocate();


}

/// This function initializes the paths depending on how they are
/// specified to be initialized in the input file.  Currently, the
/// options are :
/// CUBIC:
/// BCC:
/// FIXED:
/// FILE:
/// LEVIFLIGHT
void PathClass::InitPaths (IOSectionClass &in)
{
  SetMode (NEWMODE);
  assert(in.OpenSection ("Particles"));
  int numSpecies = NumSpecies();
  // Now initialize the Path
  for (int speciesIndex=0; speciesIndex<numSpecies; speciesIndex++)
  {
    SpeciesClass &species = *SpeciesArray(speciesIndex);
    assert(in.OpenSection("Species", speciesIndex));
    string InitPaths;
    in.ReadVar ("InitPaths", InitPaths);
    string Replicate;
    in.ReadVar ("Replicate", Replicate);
    if (InitPaths == "RANDOM") {
      cerr << "Don't know how to do RANDOM yet.\n";
      exit(1);
    }
    else if (InitPaths == "CUBIC") {
      int num = species.NumParticles;
      bool isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
      if (!isCubic) {
	cerr << "A cubic box is current required for cubic initilization\n";
	abort();
      }
      int numPerDim = (int) ceil (pow((double)num, 1.0/3.0)-1.0e-6);
      double delta = Box[0] / numPerDim;
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
	int ix, iy, iz;
	int ip = ptcl-species.FirstPtcl;
	ix = ip/(numPerDim*numPerDim);
	iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
	iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
	dVec r;
	r[0] = ix*delta;
	r[1] = iy*delta;
	r[2] = iz*delta;
	for (int slice=0; slice<NumTimeSlices(); slice++) 
	  Path(slice,ptcl) = r;
      }
    }
    else if (InitPaths == "BCC") {
      int num = species.NumParticles;
      bool isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
      if (!isCubic) {
	cerr << "A cubic box is current required for cubic initilization\n";
	abort();
      }
      int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/3.0)-1.0e-6);
      double delta = Box[0] / numPerDim;
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
	int ip = (ptcl-species.FirstPtcl)/2;
	int ix, iy, iz;
	ix = ip/(numPerDim*numPerDim);
	iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
	iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
	dVec r;
	r[0] = ix*delta-0.5*Box[0];
	r[1] = iy*delta-0.5*Box[1];
	r[2] = iz*delta-0.5*Box[2];
	if (ptcl % 2) 
	  r += 0.5*delta;
	fprintf (stderr, "BCC ptcl %d position = [%8.4f %8.4f %8.4f]\n",
		 r[0], r[1], r[2]);
	for (int slice=0; slice<NumTimeSlices(); slice++) 
	  Path(slice,ptcl) = r;
      }
    }
    else if (InitPaths == "FIXED") {
      Array<double,2> Positions;
      assert (in.ReadVar ("Positions", Positions));
      
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      for (int ptcl=species.FirstPtcl; 
	   ptcl<=species.LastPtcl; ptcl++){
	for (int slice=0; slice<NumTimeSlices(); slice++) {
	  dVec pos;
	  pos = 0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
	  Path(slice,ptcl) = pos;
	}      
      }
    }    
    else if (InitPaths == "FILE"){
      bool replicate = false;
      if (Replicate == "ON"){
        replicate = true;
        cerr << "Replicate 'ON'; Using time slice 0 for all time slices." << endl;
      }
      string pathFile;
      assert(in.ReadVar("File",pathFile));
      ReadOld(pathFile,replicate);
    }
    else if (InitPaths == "LEVIFLIGHT") {
      Array<double,2> Positions;
      assert (in.ReadVar ("Positions", Positions));
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      Array<dVec,1> R0(species.NumParticles);
      for (int ptcl=0l; ptcl<species.NumParticles; ptcl++) 
	for (int dim=0; dim<NDIM; dim++)
	  R0(ptcl)[dim] = Positions(ptcl,dim);
      NodeAvoidingLeviFlight (speciesIndex,R0);
    }
//     else if (InitPaths == "LEVIFLIGHT") {
//       int myStart, myEnd;
//       SliceRange (Communicator.MyProc(), myStart, myEnd);
//       Array<double,2> Positions;
//       assert (in.ReadVar ("Positions", Positions));
//       assert (Positions.rows() == species.NumParticles);
//       assert (Positions.cols() == species.NumDim);
//       Array<dVec,1> flight(TotalNumSlices+1);
//       for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
// 	for (int dim=0; dim<NDIM; dim++) {
// 	  flight(0)[dim] = Positions(ptcl-species.FirstPtcl,dim);
// 	  flight(TotalNumSlices)[dim] = Positions(ptcl-species.FirstPtcl,dim);
// 	}
// 	LeviFlight (flight, species.lambda);
// // 	for (int i=0; i<MyNumSlices; i++)
// // 	  Path(i, ptcl) = flight(i-RefSlice);
// 	for (int slice=myStart; slice<=myEnd; slice++)
// 	  Path(slice-myStart, ptcl) = flight(slice);
//       }
//     }

    else {
      cerr << "Unrecognize initialization strategy " 
	   << InitPaths << endl;
      abort();
    }
    in.CloseSection(); // Species
//   jgadd: get numMol
    if (speciesIndex == 0){  
      numMol = species.NumParticles;
    }
  }
  in.CloseSection(); // "Particles"

// jgadd: correct entries where necessary
//  cerr << "PRINT MolRef corrected" << endl;
  for(int m = 0;m < MolRef.size(); m++){
    int ref = MolRef(m);
    while(ref >= numMol){
      ref -= numMol;
    MolRef(m) = ref;
    }
    cerr << m << " " << MolRef(m) << endl;
  } 
  string openSpeciesName;
  if (OpenPaths) {
    assert(in.ReadVar("OpenSpecies",openSpeciesName));
    OpenSpeciesNum=SpeciesNum(openSpeciesName);
    InitOpenPaths();
  }  

  //Everything needs to be accepted
  Path.AcceptCopy();
  Permutation.AcceptCopy();
  Rho_k.AcceptCopy();
  BroadcastRefPath();
  RefPath.AcceptCopy();
  Weight.AcceptCopy();
}


inline bool Include(dVec k)
{
  //  assert (NDIM == 3);
  if (k[0] > 0.0)
    return true;
  else if ((k[0]==0.0) && (k[1]>0.0))
    return true;
  else if ((NDIM==3) && ((k[0]==0.0) && (k[1]==0.0) && (k[2] > 0.0)))
    return true;
  else
    return false;
}
 

void PathClass::Allocate()
{

  assert(TotalNumSlices>0);
  int myProc=Communicator.MyProc();
  int numProcs=Communicator.NumProcs();
  /// Everybody gets the same number of time slices if possible.
  /// Otherwise the earlier processors get the extra one slice  
  /// until we run out of extra slices. The last slice on processor i
  /// is the first slice on processor i+1.
  MyNumSlices=TotalNumSlices/numProcs+1+(myProc<(TotalNumSlices % numProcs));
  
  // Initialize reference slice position to be 0 on processor 0.  
  RefSlice = 0;
    
  int numParticles = 0;

  /// Set the particle range for the new species
  for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
    SpeciesArray(speciesNum)->FirstPtcl = numParticles;
    numParticles=numParticles + SpeciesArray(speciesNum)->NumParticles;
    SpeciesArray(speciesNum)->LastPtcl= numParticles-1;
  }
  Path.resize(MyNumSlices,numParticles+OpenPaths);
  RefPath.resize(numParticles);
  Permutation.resize(numParticles+OpenPaths);
  SpeciesNumber.resize(numParticles+OpenPaths);
  DoPtcl.resize(numParticles+OpenPaths);
  /// Assign the species number to the SpeciesNumber array
  for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
    for (int i=SpeciesArray(speciesNum)->FirstPtcl; 
	 i<= SpeciesArray(speciesNum)->LastPtcl; i++)
      SpeciesNumber(i) = speciesNum;
  }
  //Sets to the identity permutaiton 
  for (int ptcl=0;ptcl<Permutation.size();ptcl++){
    Permutation(ptcl) = ptcl;
  }
  if (LongRange) {
#if NDIM==3    
    SetupkVecs3D();
#endif
#if NDIM==2
    SetupkVecs2D();
#endif
    Rho_k.resize(MyNumSlices, NumSpecies(), kVecs.size());
  }
  // jgadd
  MolRef.resize(numParticles);
//  cerr << "PRINT MolRef Initialized" << endl;
  for(int q = 0;q < MolRef.size();q++){
    MolRef(q) = q;
//    cerr << q << " " << MolRef(q) << endl;
  }
}

void PathClass::SetupkVecs2D()
{

  for (int i=0; i<NDIM; i++)
    kBox[i] = 2.0*M_PI/Box[i];
  
  int numVecs=0;

  assert (NDIM == 2);

  for (int i=0;i<NDIM;i++){
    MaxkIndex[i]= (int) ceil(1.1*kCutoff/kBox[i]);
    // cerr << "MaxkIndex[" << i << "] = " << MaxkIndex[i] << endl;
  }


  dVec k;
  TinyVector<int,NDIM> ki;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      //      for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
      //	k[2] = iz*kBox[2];
      if ((dot(k,k)<kCutoff*kCutoff) && Include(k))
	numVecs++;
      //}
    }
  }
  kIndices.resize(numVecs);
  cerr << "kCutoff = " << kCutoff << endl;
  cerr << "Number of kVecs = " << numVecs << endl;
  kVecs.resize(numVecs);
  for (int i=0; i<NDIM; i++)
    C[i].resize(2*MaxkIndex[i]+1);
  numVecs = 0;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    ki[0]= ix+MaxkIndex[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      ki[1]= iy+MaxkIndex[1];
      //      for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
      //	k[2] = iz*kBox[2];
      //	ki[2]= iz+MaxkIndex[2];
	if ((dot(k,k)<kCutoff*kCutoff) && Include(k)) {
	  cerr<<"This k vec is "<<k[0]<<" "<<k[1]<<endl;
	  kVecs(numVecs) = k;
	  kIndices(numVecs)=ki;
	  numVecs++;
	}
	//      }
    }
  }
  SortRhoK();
}


///Puts in the vector MagKInt a number that corresponds to the sorted
///order of teh magnitude of the k vectors. Doesn't actually change
///what's in kVec. Currently only used for the reading of David's long
///range class. 
void PathClass::SortRhoK()
{
  /////  Array<double,1> MagK;
  //  Array<int,1> MagKint;
  MagK.resize(kVecs.size());
  MagKint.resize(kVecs.size());
  for (int kVec=0;kVec<kVecs.size();kVec++){
    MagK(kVec)=sqrt(kVecs(kVec)[0]*kVecs(kVec)[0]+
			kVecs(kVec)[1]*kVecs(kVec)[1]);
  }
  int smallIndex=0;
 
  for (int counter=0;counter<MagKint.size();counter++){
    MagKint(counter)=-1;
  }



  for (int currentNum=0;currentNum<kVecs.size();currentNum++){
    for (int counter=0;counter<MagK.size();counter++){
      if (MagKint(counter)==-1 && MagK(counter)<MagK(smallIndex)){
	smallIndex=counter;
      }
    }
    for (int counter=0;counter<MagK.size();counter++){
      if (abs(MagK(smallIndex)-MagK(counter))<1e-4){
	MagKint(counter)=currentNum;
      }
    }
    smallIndex=-1;
    for (int counter=0;counter<MagK.size();counter++){
      if (MagKint(counter)==-1)
	smallIndex=counter;
    }
    if (smallIndex==-1)
      currentNum=kVecs.size()+1;
  }
  for (int counter=0;counter<MagKint.size();counter++){
    cerr<<"My mag K int is "<<MagKint(counter)<<endl;
  }

}


void PathClass::SetupkVecs3D()
{

  for (int i=0; i<NDIM; i++)
    kBox[i] = 2.0*M_PI/Box[i];
  
  int numVecs=0;

  assert (NDIM == 3);

  for (int i=0;i<NDIM;i++){
    MaxkIndex[i]= (int) ceil(1.1*kCutoff/kBox[i]);
    // cerr << "MaxkIndex[" << i << "] = " << MaxkIndex[i] << endl;
  }


  dVec k;
  TinyVector<int,NDIM> ki;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
	k[2] = iz*kBox[2];
	if ((dot(k,k)<kCutoff*kCutoff) && Include(k))
	  numVecs++;
      }
    }
  }
  kIndices.resize(numVecs);
  cerr << "kCutoff = " << kCutoff << endl;
  cerr << "Number of kVecs = " << numVecs << endl;
  kVecs.resize(numVecs);
  for (int i=0; i<NDIM; i++)
    C[i].resize(2*MaxkIndex[i]+1);
  numVecs = 0;
  for (int ix=-MaxkIndex[0]; ix<=MaxkIndex[0]; ix++) {
    k[0] = ix*kBox[0];
    ki[0]= ix+MaxkIndex[0];
    for (int iy=-MaxkIndex[1]; iy<=MaxkIndex[1]; iy++) {
      k[1] = iy*kBox[1];
      ki[1]= iy+MaxkIndex[1];
      for (int iz=-MaxkIndex[2]; iz<=MaxkIndex[2]; iz++) {
	k[2] = iz*kBox[2];
	ki[2]= iz+MaxkIndex[2];
	if ((dot(k,k)<kCutoff*kCutoff) && Include(k)) {
	  kVecs(numVecs) = k;
	  kIndices(numVecs)=ki;
	  numVecs++;
	}
      }
    }
  }
}


void PathClass::CalcRho_ks_Slow(int slice, int species)
{
  for (int ki=0; ki<kVecs.size(); ki++) {
    complex<double> rho;
    rho = 0.0;
    for (int ptcl=Species(species).FirstPtcl; 
	 ptcl <= Species(species).LastPtcl; ptcl++) {
      const dVec &r = (*this)(slice, ptcl);
      double phase = dot(r, kVecs(ki));
      rho += complex<double> (cos(phase), sin(phase));
    }
    Rho_k(slice, species, ki) = rho;
  }
}

void PathClass::CalcRho_ks_Fast(int slice,int species)
{
  //  cerr<<"Beginning the calcrhok stuff"<<endl;
  // Zero out Rho_k array
  for (int ki=0;ki<kIndices.size();ki++)
    Rho_k(slice,species,ki)=0.0;

  for (int ptcl=Species(species).FirstPtcl; 
       ptcl <= Species(species).LastPtcl; ptcl++) {
    // First, compute C arrays
    const dVec &r = (*this)(slice, ptcl);
    for (int dim=0;dim<NDIM;dim++){
      complex<double> tempC;
      double phi = r[dim] * kBox[dim];
      tempC=complex<double>(cos(phi), sin(phi));
      C[dim](MaxkIndex[dim]) = 1.0;
      for (int n=1; n<=MaxkIndex[dim]; n++) {
	C[dim](MaxkIndex[dim]+n) = tempC * C[dim](MaxkIndex[dim]+n-1);
	C[dim](MaxkIndex[dim]-n) = conj(C[dim](MaxkIndex[dim]+n));
      }
    }
    // Now, loop over k-vector indices;
    for (int ki=0; ki<kIndices.size(); ki++) {
      const TinyVector<int,NDIM> &kIndex = kIndices(ki);
      //#ifdef THREE_D
#if NDIM==3
      Rho_k(slice,species,ki) += 
	C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
      //#ifdef TWO_D
#if NDIM==2
      Rho_k(slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1]);
#endif
    }
  }
  //  cerr<<"ending the calcrhok stuff"<<endl;
}
      


void PathClass::MoveJoin(int oldJoin, int newJoin)
{
  //  cerr<<"My time slices is "<<NumTimeSlices()<<endl;
  //  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    //    cerr<<Permutation(ptcl)<<endl;
  //  }
  //  cerr<<oldJoin<<" "<<newJoin<<" "<<endl;
  //  cerr<<"Starting"<<endl;
  //  cerr<<Path(OpenLink,OpenPtcl)<<endl;
  bool swappedAlready=false;
  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	Path[OLDMODE](timeSlice,ptcl)=Path[NEWMODE](timeSlice,Permutation(ptcl));
	if (timeSlice==(int)OpenLink && Permutation(ptcl)==(int)OpenPtcl && !swappedAlready){
	  OpenPtcl[OLDMODE]=ptcl;
	  OpenPtcl[NEWMODE]=ptcl;
	  swappedAlready=true;
	}
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){ 
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	Path[NEWMODE](timeSlice,ptcl)=Path[OLDMODE](timeSlice,ptcl);

      }
    }
  }
  else if (oldJoin>newJoin){
    //  else if (oldJoin>=newJoin){//CHANGED!
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	if (timeSlice==(int)OpenLink && ptcl==(int)OpenPtcl && !swappedAlready){
	  OpenPtcl[OLDMODE]=Permutation(ptcl);
	  OpenPtcl[NEWMODE]=Permutation(ptcl);
	  swappedAlready=true;
	}
	Path[OLDMODE](timeSlice,Permutation(ptcl))=Path[NEWMODE](timeSlice,ptcl);
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	Path[NEWMODE](timeSlice,ptcl)=Path[OLDMODE](timeSlice,ptcl);
	
      }
    }
  }
  //  cerr<<Path(OpenLink,OpenPtcl)<<endl;
  //  cerr<<"Ending"<<endl;
}



void PathClass::AcceptCopy(int startSlice,int endSlice, 
				  const Array <int,1> &activeParticles)
{
  Weight.AcceptCopy();
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[OLDMODE](Range(startSlice, endSlice), ptcl) = 
      Path[NEWMODE](Range(startSlice, endSlice), ptcl);
    Permutation.AcceptCopy(ptcl);
  }

  if (OpenPaths){
    OpenPtcl.AcceptCopy();
    OpenLink.AcceptCopy();
    for (int counter=0;counter<NumTimeSlices();counter++){
      Path[OLDMODE](counter,NumParticles())=Path[NEWMODE](counter,NumParticles());
    }
    //    Path[OLDMODE](Range(startSlice,endSlice),NumParticles())=
    //      Path[NEWMODE](Range(startSlice,endSlice),NumParticles());
  }
    
  
}

void PathClass::RejectCopy(int startSlice,int endSlice, 
				  const Array <int,1> &activeParticles)
{
  Weight.RejectCopy();
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[NEWMODE](Range(startSlice, endSlice), ptcl) = 
      Path[OLDMODE](Range(startSlice, endSlice), ptcl);
    Permutation.RejectCopy(ptcl);
  }

  if (OpenPaths){
    OpenPtcl.RejectCopy();
    OpenLink.RejectCopy();
    Path[NEWMODE](Range(startSlice,endSlice),NumParticles())=
      Path[OLDMODE](Range(startSlice,endSlice),NumParticles());
  }

}


void PathClass::ShiftData(int slicesToShift)
{

  ShiftPathData(slicesToShift);
  if (LongRange)
    ShiftRho_kData(slicesToShift);
  OpenLink.AcceptCopy(); ///the open link has changed and you want to accept it
  RefSlice += slicesToShift;
  while (RefSlice >= TotalNumSlices)
    RefSlice -= TotalNumSlices;
  while (RefSlice < 0)
    RefSlice += TotalNumSlices;
  //  cerr<<"My ref slice at particle 0 is "<<Path(RefSlice,0)<<endl;
  if (OpenPaths){
    cerr<<"Here my open link is "<<OpenLink<<endl;
    int openLinkOld=(int)OpenLink;
    OpenLink=RefSlice;
    if ((int)OpenLink==0){
      OpenLink=NumTimeSlices()-1;
    }
    //  cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    ////    cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<" "<<slicesToShift<<endl;
    ////    if (openLinkOld!=(int)OpenLink)
    ////      cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<" "<<slicesToShift<<endl;
    ////    assert(openLinkOld==(int)OpenLink);
    OpenLink[OLDMODE]=OpenLink[NEWMODE];
  }
}

void PathClass::ShiftRho_kData(int slicesToShift)
{

  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numSlices  = Rho_k.extent(0);
  int numSpecies = Rho_k.extent(1);
  int numk       = Rho_k.extent(2);
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  /// First shifts the data in the A copy left or right by the
  /// appropriate amount 
  if (slicesToShift>0) {
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int species=0; species<numSpecies; species++)
	for (int ki=0; ki<numk; ki++)
	  Rho_k[0](slice,species,ki)=Rho_k[0](slice-slicesToShift,species,ki);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int species=0; species<numSpecies; species++)
	for (int ki=0; ki<numk; ki++)
	  Rho_k[0](slice,species,ki)=Rho_k[0](slice-slicesToShift,species,ki);
  }
  

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numSpecies*numk;
  Array<complex<double>,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0) {
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++)
      for (int species=0; species<numSpecies; species++)
	for (int ki=0; ki<numk; ki++) {
	  /// If shifting forward, don't send the last time slice (so always)
	  /// send slice-1
	  sendBuffer(buffIndex)=Rho_k[1](slice-1,species,ki);
	  buffIndex++;
	}
  }
  else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++)
      for (int species=0; species<numSpecies; species++)
	for (int ki=0; ki<numk; ki++) {
	  /// If shifting backward, don't send the first time slice (so always)
	  /// send slice+1
	  sendBuffer(buffIndex)=Rho_k[1](slice+1,species,ki);
	  buffIndex++;
	}
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
  if (slicesToShift>0)
    startSlice=0;
  else 
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++)
    for (int species=0; species<numSpecies; species++)
      for (int ki=0; ki<numk; ki++) {
	Rho_k[0](slice,species,ki)=receiveBuffer(buffIndex);
	buffIndex++;
      }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int species=0; species<numSpecies; species++)
      for (int ki=0; ki<numk; ki++) 
	Rho_k[1](slice,species,ki) = Rho_k[0](slice,species,ki);

  // And we're done! 
      
}

void PathClass::ShiftPathData(int slicesToShift)
{
  //  cerr<<"Slices to shift are "<<slicesToShift<<endl;
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numPtcls=NumParticles()+OpenPaths;
  int numSlices=NumTimeSlices();
  assert(abs(slicesToShift)<numSlices);
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  ///First shifts the data in the A copy left 
  ///or right by the appropriate amount   
  if (slicesToShift>0){
    for (int slice=numSlices-1; slice>=slicesToShift;slice--)
      for (int ptcl=0;ptcl<numPtcls;ptcl++)
	Path[NEWMODE](slice,ptcl) = Path[NEWMODE](slice-slicesToShift,ptcl);
  }
  else {
    for (int slice=0; slice<numSlices+slicesToShift;slice++)
      for (int ptcl=0;ptcl<numPtcls;ptcl++)
	Path[NEWMODE](slice,ptcl) = Path[NEWMODE](slice-slicesToShift,ptcl);
  }
  

  /// Now bundle up the data to send to adjacent processor
  int bufferSize=abs(slicesToShift)*numPtcls;
  Array<dVec,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startSlice;
  int buffIndex=0;
  if (slicesToShift>0){
    startSlice=numSlices-slicesToShift;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int ptcl=0;ptcl<numPtcls;ptcl++){
	///If shifting forward, don't send the last time slice (so always)
	///send slice-1
	sendBuffer(buffIndex)=Path[OLDMODE](slice-1,ptcl);
	buffIndex++;
      }
    }
  }
  else {
    startSlice=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      for (int ptcl=0;ptcl<numPtcls;ptcl++){
	///If shifting backward, don't send the first time slice (so always)
	///send slice+1
	sendBuffer(buffIndex)=Path[OLDMODE](slice+1,ptcl);
	buffIndex++;
      }
    }
  }

  /// Send and receive data to/from neighbors.
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
  if (slicesToShift>0)
    startSlice=0;
  else 
    startSlice=numSlices+slicesToShift;

  /// Copy the data into the A copy
  buffIndex=0;
  for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
    for (int ptcl=0;ptcl<numPtcls;ptcl++){
      Path[NEWMODE](slice,ptcl)=receiveBuffer(buffIndex);
      buffIndex++;
    }
  }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int slice=0; slice<numSlices; slice++)
    for (int ptcl=0; ptcl<numPtcls; ptcl++)
      Path[OLDMODE](slice,ptcl) = Path[NEWMODE](slice,ptcl);

  // And we're done!
  if (OpenPaths){ //only works for serialo positive slices
    if ((int)OpenLink+slicesToShift >= NumTimeSlices())
      OpenLink=OpenLink+1;
    if ((int)OpenLink+slicesToShift <0)
      OpenLink=OpenLink-1;
    OpenLink=((int)OpenLink+slicesToShift+NumTimeSlices()) % NumTimeSlices();
    if ((int)OpenLink==0){
      cerr<<"equal to 0"<<endl;
      OpenLink=NumTimeSlices()-1;
    }

    //    OpenLink=RefSlice;
    //    if ((int)OpenLink==0){
    //      OpenLink=NumTimeSlices()-1;
    //    }
    //    cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    //    if (openLinkOld!=(int)OpenLink)
    //      cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    //    assert(openLinkOld==(int)OpenLink);
    OpenLink[OLDMODE]=OpenLink[NEWMODE];
  }
  //  cerr<<Path(OpenLink,NumParticles())<<endl;
}


/// This function sends the reference slice from whatever processor is
/// currently holding it to all the other processors.  Note that only
/// the NEW copy is updated.  If we want the OLD copy to be updated,
/// we must use RefPath.AcceptCopy().
void PathClass::BroadcastRefPath()
{
  SetMode (NEWMODE);
  // Figure which processor has the reference slice
  int procWithRefSlice = SliceOwner (RefSlice);
  
  Array<dVec,1> buffer(NumParticles());
  if (procWithRefSlice == Communicator.MyProc()) {
    int myStart, myEnd;
    SliceRange (Communicator.MyProc(), myStart, myEnd);
    for (int ptcl=0; ptcl<NumParticles(); ptcl++)
      buffer(ptcl) = Path(RefSlice-myStart, ptcl);
  }
  
  //  cerr << "procWithRefSlice = " << procWithRefSlice << endl;

  // Do broadcast
  Communicator.Broadcast(procWithRefSlice, buffer);

  // Now, all processors have reference slice.  Note that the 
  // labeling of particles may not be right since with haven't
  // propogated the permuation matrices.
  
  // Now, copy into path
  for (int ptcl=0; ptcl<NumParticles(); ptcl++)
    RefPath(ptcl) = buffer(ptcl);
}


void PathClass::TotalPermutation(Array<int,1> &permVec)
{
  int myProc = Communicator.MyProc();
  int numProcs = Communicator.NumProcs();
  int numPtcls = NumParticles();
  permVec.resize(Permutation.size());
  for (int i=0; i<permVec.size(); i++)
    permVec(i) = Permutation(i);
  Array<int,2> permMat(numProcs, permVec.size());

  // First, collect all the individual permutation vectors
  Communicator.Gather (permVec, permMat, 0);
  // Now apply them in sequence.
  for (int pi=0; pi < numPtcls; pi++) {
    int ptcl = pi;
    for (int proc=0; proc<numProcs; proc++)
      ptcl = permMat(proc, ptcl);
    permVec(pi) = ptcl;
  }
}
