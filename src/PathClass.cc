#include "PathClass.h"

void PathClass::Read (IOSectionClass &inSection)
{
  SetMode (NEWMODE);
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
	ix = ptcl/(numPerDim*numPerDim);
	iy = (ptcl-(ix*numPerDim*numPerDim))/numPerDim;
	iz = ptcl - ix*numPerDim*numPerDim - iy*numPerDim;
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
	int ip = ptcl/2;
	int ix, iy, iz;
	ix = ip/(numPerDim*numPerDim);
	iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
	iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
	dVec r;
	r[0] = ix*delta;
	r[1] = iy*delta;
	r[2] = iz*delta;
	if (ptcl % 2) 
	  r += 0.5*delta;
	for (int slice=0; slice<NumTimeSlices(); slice++) 
	  Path(slice,ptcl) = r;
      }
    }
    else if (InitPaths == "FIXED") {
      Array<double,2> Positions;
      assert (inSection.ReadVar ("Positions", Positions));
      
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      for (int ptcl=species.FirstPtcl; 
	   ptcl<=species.LastPtcl; ptcl++){
	for (int slice=0; slice<NumTimeSlices(); slice++) {
	  cerr<<ptcl;
	  dVec pos;
	  pos = 0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
	  Path(slice,ptcl) = pos;
	}      
      }
    }
    else {
      cerr << "Unrecognize initialization strategy " 
	   << InitPaths << endl;
    }
    inSection.CloseSection();
  }
  inSection.CloseSection(); // "Particles"
  //Everything needs to be accepted
  Path.AcceptCopy();
  Permutation.AcceptCopy();
  Rho_k.AcceptCopy();
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
//   cerr<<"Numprocs is "<<numProcs<<endl;
//   cerr<<"mynumslices: "<<MyNumSlices<<endl;
  
  
  int numParticles = 0;

  /// Set the particle range for the new species
  for (int speciesNum=0;speciesNum<SpeciesArray.size();speciesNum++){
    SpeciesArray(speciesNum)->FirstPtcl = numParticles;
    numParticles=numParticles + SpeciesArray(speciesNum)->NumParticles;
    SpeciesArray(speciesNum)->LastPtcl= numParticles-1;
  }
  cerr<<"my number of particles is "<<numParticles<<endl;
  Path.resize(MyNumSlices,numParticles);
  Permutation.resize(numParticles);
  SpeciesNumber.resize(numParticles);
  DoPtcl.resize(numParticles);
  cerr<<"Permutation size is "<<Permutation.size()<<endl;
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
  if (LongRange){
    SetupkVecs();
    Rho_k.resize(MyNumSlices, NumSpecies(),kVecs.size());
  }
  
}


void PathClass::SetupkVecs()
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
#ifdef THREE_D
      Rho_k(slice,species,ki) += 
	C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
#ifdef TWO_D
      Rho_k(slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1]);
#endif
    }
  }
}
      


void PathClass::MoveJoin(int oldJoin, int newJoin)
{
  //  cerr<<"My time slices is "<<NumTimeSlices()<<endl;
  //  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    //    cerr<<Permutation(ptcl)<<endl;
  //  }

  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	Path[OLDMODE](timeSlice,ptcl)=Path[NEWMODE](timeSlice,Permutation(ptcl));
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
}



void PathClass::AcceptCopy(int startSlice,int endSlice, 
				  const Array <int,1> &activeParticles)
{

  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[OLDMODE](Range(startSlice, endSlice), ptcl) = 
      Path[NEWMODE](Range(startSlice, endSlice), ptcl);
    Permutation.AcceptCopy(ptcl);
  }
}

void PathClass::RejectCopy(int startSlice,int endSlice, 
				  const Array <int,1> &activeParticles)
{

  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[NEWMODE](Range(startSlice, endSlice), ptcl) = 
      Path[OLDMODE](Range(startSlice, endSlice), ptcl);
    Permutation.RejectCopy(ptcl);
  }
}


void PathClass::ShiftData(int slicesToShift)
{
  ShiftPathData(slicesToShift);
  if (LongRange)
    ShiftRho_kData(slicesToShift);
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
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc, sendProc;
  int numPtcls=NumParticles();
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
}
