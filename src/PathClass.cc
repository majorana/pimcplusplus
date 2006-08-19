/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "PathClass.h"
#include "Actions/ActionsClass.h"

double PathClass::MinImageDistance(dVec v1, dVec v2)
{
  dVec disp = v2-v1;
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  double dist = sqrt(dot(disp,disp));
  return dist;
}

dVec PathClass::MinImageDisp(dVec v1, dVec v2)
{
  dVec disp = v2-v1;
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  return disp;
}


const dVec& 
PathClass::GetOpenTail()
{
  return Path(OpenLink,NumParticles());
}

const dVec
PathClass::ReturnOpenHead()
{
  return Path(OpenLink,OpenPtcl);
}

// //NOTE: Be careful.  
// const dVec& 
// PathClass::GetOpenHead()
// {
//   return Path(OpenLink,OpenPtcl);
// }

void 
PathClass::SetHead(const dVec &r,int join)
{
  (*this)(OpenLink, OpenPtcl) = r;
  if (OpenLink==NumTimeSlices()-1){
    if (join==NumTimeSlices()-1)
      (*this)(0,Permutation(OpenPtcl))=r;
    else 
      (*this)(0,OpenPtcl)=r;
  }
}

void 
PathClass::SetTail(const dVec &r)
{
  (*this)(OpenLink, NumParticles()) = r;
}




void PathClass::RefDistDisp (int slice, int refPtcl, int ptcl,
			     double &dist, dVec &disp)
{
  disp = Path(slice, ptcl)- RefPath(refPtcl);
  
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
    if (!(-Box(i)/2.0<=disp(i))){
      perr<<"ERROR: "<<Box(i)<<" "<<disp(i)<<" "
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
      perr<<DBdisp(i)<<" "<<disp(i)<<endl;
    }
    //    assert (fabs(DBdisp(i)-disp(i)) < 1.0e-12);
  }
#endif
}



void
PathClass::SetIonConfig(int config)
{
  ConfigNum = config;
  int first = Species(IonSpecies).FirstPtcl;
  
  SetMode(OLDMODE);
  for (int ptcl=0; ptcl<IonConfigs[config].size(); ptcl++) 
    for (int slice=0; slice<NumTimeSlices(); slice++)
      SetPos(slice, ptcl+first, IonConfigs[config](ptcl));

  SetMode(NEWMODE);
  for (int ptcl=0; ptcl<IonConfigs[config].size(); ptcl++) 
    for (int slice=0; slice<NumTimeSlices(); slice++)
      SetPos(slice, ptcl+first, IonConfigs[config](ptcl));
}

void PathClass::Read (IOSectionClass &inSection)
{
  inSection.ReadVar("FunnyCoupling",FunnyCoupling);
  SetMode(OLDMODE);
  NowOpen=false;
  HeadSlice=0;
  SetMode (NEWMODE);
  NowOpen=false;
  HeadSlice=0;
  Weight=1;
  double tempExistsCoupling;
  if (!inSection.ReadVar("ExistsCoupling",tempExistsCoupling)){
    ExistsCoupling=-1.0;
  }
  else{
    ExistsCoupling=tempExistsCoupling;
  }
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
    periodic[i] = tempPeriodic(i);
  }
  SetPeriodic (periodic);
  if (needBox) {
    assert(inSection.ReadVar ("Box", tempBox));
    perr << "Using periodic boundary conditions.\n";
    assert(tempBox.size()==NDIM);
    //HACK!
    // //    int desiredNumParticles;
    // //    assert(inSection.ReadVar("NumParticles",desiredNumParticles));
    // //    double desiredBoxVol=(tempBox(0)*tempBox(1)*tempBox(2));
    // //    double desiredDensity=
    // //      desiredNumParticles/desiredBoxVol;
    // //    ScaleBox=
    // //      pow((double)(desiredNumParticles-MyClone)/(double)(desiredNumParticles),1.0/3.0);
    //END HACK
    ScaleBox=1.0;

    for (int counter=0;counter<tempBox.size();counter++)
      Box(counter)=tempBox(counter)*ScaleBox;
    SetBox (Box);
  }
  else 
    perr << "Using free boundary conditions.\n";
  OrderN=false;
  inSection.ReadVar("OrderN",OrderN);
  if (!inSection.ReadVar("OpenLoops",OpenPaths))
    OpenPaths=false;

#ifndef OPEN_LOOPS
  if (OpenPaths) {
    cerr << "OpenPaths are not enabled in the code!\n"
	 << "Reconfigure with --enable-open and \n" 
	 << "  make clean; make ...\nAborting.\n";
    abort();
  }
#endif

  // Read in the k-space radius.  If we don't have that,
  // we're not long-ranged.
  LongRange = inSection.ReadVar("kCutoff", kCutoff);
  if (!(inSection.ReadVar("DavidLongRange",DavidLongRange))){
    DavidLongRange=false;
  }
  if (DavidLongRange)
    perr<<"I am doing DAVID LONG RANGE!"<<endl;
  assert(inSection.OpenSection("Particles"));
  int numSpecies = inSection.CountSections ("Species");
  perr<<"we have this many sections: "<<numSpecies<<endl;
   // First loop over species and read info about species
  for (int Species=0; Species < numSpecies; Species++) {
    inSection.OpenSection("Species", Species);
    SpeciesClass *newSpecies = ReadSpecies (inSection);
		doMol = newSpecies->AssignMoleculeIndex;
    inSection.CloseSection(); // "Species"
    bool manyParticles=false;
    inSection.ReadVar("ManyParticles",manyParticles);
    if (manyParticles){
      //UGLY HACK!
      newSpecies->NumParticles=newSpecies->NumParticles-MyClone;
      //UGLY HACK!
    }
    AddSpecies (newSpecies);
  }
  inSection.CloseSection(); // Particles
  // Now actually allocate the path
  Allocate();


  /// Read to see if we are using correlated sampling
  string ionSpecies;
  if (inSection.ReadVar ("IonSpecies", ionSpecies)) {
    CorrelatedSampling = true;
    IonSpecies = SpeciesNum(ionSpecies);
    assert(IonSpecies!=-1);
    Array<double,2> ionConfigA, ionConfigB;
    assert(inSection.ReadVar("IonConfigA", ionConfigA));
    assert(inSection.ReadVar("IonConfigB", ionConfigB));
    assert(ionConfigA.extent(0) == ionConfigB.extent(0));
    assert(ionConfigA.extent(1) == NDIM);
    assert(ionConfigB.extent(1) == NDIM);
    IonConfigs[0].resize(ionConfigA.extent(0));
    IonConfigs[1].resize(ionConfigB.extent(0));
    for (int i=0; i<ionConfigA.extent(0); i++)
      for (int j=0; j<3; j++) {
	IonConfigs[0](i)[j] = ionConfigA(i,j);
	IonConfigs[1](i)[j] = ionConfigB(i,j);
      }
    cerr << "IonConfigs[0] = " << IonConfigs[0] << endl;
    cerr << "IonConfigs[1] = " << IonConfigs[1] << endl;
    SetIonConfig(0);
  }
  NowOpen=true;

  /// Checking rounding mode
  bool roundOkay = true;
  for (int i=0; i<1000; i++) {
    double x = 10.0*drand48()-5.0;
    if (nearbyint(x) != round(x))
      roundOkay = false;
  }
  if (!roundOkay) {
    cerr << "Rounding mode is not set to ""round"".  Aborting!\n";
    abort();
  }

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
    cerr << "  species " << (speciesNum+1) << " of " 
	 << SpeciesArray.size() << endl;
    int N = SpeciesArray(speciesNum)->NumParticles;
    cerr << "		I have N = " << N << " and doMol is " << doMol << endl;
    int first = numParticles;
    SpeciesArray(speciesNum)->FirstPtcl = first;
    numParticles=numParticles + N;
    SpeciesArray(speciesNum)->LastPtcl = first + N-1;;
    SpeciesArray(speciesNum)->Ptcls.resize(N);
    
    int foundAt;
    int prevIndex;
    if(doMol){
      string newMol = SpeciesArray(speciesNum)->molecule;
      cerr << "		Looking for string " << newMol << endl;
      bool foundMol = false;
      for(int m=0; m<MoleculeName.size(); m++){
				if(newMol == MoleculeName[m]){
				  foundMol = true;
				  assert(N/SpeciesArray(speciesNum)->formula == MoleculeNumber[m]);
				  foundAt = m;
				  cerr << "		Found at " << foundAt 
				       << " of " << MoleculeName.size() << endl;
				}
      }
      if(!foundMol){
				MoleculeName.push_back(newMol);
				MoleculeNumber.push_back(N/SpeciesArray(speciesNum)->formula);
				int sum = 0;
				for(int s=0; s<(MoleculeNumber.size()-1); s++)
					sum += MoleculeNumber[s];
				offset.push_back(sum);
				foundAt = MoleculeName.size() - 1;
				cerr << "		Added " << newMol << " at " << foundAt << endl;
      }
      
      prevIndex = 0;
      for(int i=0; i<foundAt; i++){
				prevIndex += MoleculeNumber[i];
      }
      //cerr << "		set prevIndex " << prevIndex << endl;
      
      //map<string,int>::iterator findName = (MoleculeMap.find(newMol));
      //if(findName != MoleculeMap.end()){
      //	assert(N/SpeciesArray(speciesNum)->formula == MoleculeNumber[m]);
      //}
      //else{
      //	MoleculeNumber.push_back(N/SpeciesArray(speciesNum)->formula);
      //	MoleculeMap[newMol] = MoleculeNumber.size();
      //}
      //cerr << "		Doing resizeAndPreserve...";
      MolRef.resizeAndPreserve(MolRef.size() + N);
      //cerr << " done" << endl;
      numMol = MoleculeNumber[0];
    }
    for (int i=0; i < N; i++){
      SpeciesArray(speciesNum)->Ptcls(i) = i+first;
      if(doMol){
	//cerr << "		" << i << " prevIndex " << prevIndex << ", foundAt " << foundAt << ", Mol.Num. " << MoleculeNumber[foundAt];
	//cerr  << " so mod is " << i%MoleculeNumber[foundAt] << endl;
				MolRef(i+first) = prevIndex + i%MoleculeNumber[foundAt]; // need to assign myMolecule
	//cerr << "Assigned MolRef " << MolRef(i+first) << " to ptcl " << i+first << endl;
      }
    }
  }
  Path.resize(MyNumSlices,numParticles+OpenPaths);
  RefPath.resize(numParticles+OpenPaths);
  Permutation.resize(numParticles+OpenPaths);
  
  /// Initilize first and last slice for all "real" particles to the
  /// first and last slices, respectively
  FirstSlice.resize(numParticles+OpenPaths);
  LastSlice.resize(numParticles+OpenPaths);
  for (int ptcl=0; ptcl<numParticles; ptcl++) {
    SetMode(OLDMODE);
    FirstSlice(ptcl) = 0;
    LastSlice(ptcl) = NumTimeSlices()-1;
    SetMode(NEWMODE);
    FirstSlice(ptcl) = 0;
    LastSlice(ptcl) = NumTimeSlices()-1;
  }
  /// Set all nonexistant particles' first and last slices to -1
  for (int ptcl=numParticles; ptcl<numParticles+OpenPaths; ptcl++) {
    SetMode(OLDMODE);
    FirstSlice(ptcl) = -1;
    LastSlice(ptcl)  = -1;
    SetMode(NEWMODE);
    FirstSlice(ptcl) = -1;
    LastSlice(ptcl)  = -1;
  }

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
}

void PathClass::SetupkVecs2D()
{

  for (int i=0; i<NDIM; i++)
    kBox[i] = 2.0*M_PI/Box[i];
  
  int numVecs=0;

  assert (NDIM == 2);

  for (int i=0;i<NDIM;i++){
    MaxkIndex[i]= (int) ceil(1.1*kCutoff/kBox[i]);
    // perr << "MaxkIndex[" << i << "] = " << MaxkIndex[i] << endl;
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
  perr << "kCutoff = " << kCutoff << endl;
  perr << "Number of kVecs = " << numVecs << endl;
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
	  perr<<"This k vec is "<<k[0]<<" "<<k[1]<<endl;
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
    perr<<"My mag K int is "<<MagKint(counter)<<endl;
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
    // perr << "MaxkIndex[" << i << "] = " << MaxkIndex[i] << endl;
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
  perr << "kCutoff = " << kCutoff << endl;
  perr << "Number of kVecs = " << numVecs << endl;
  perr << "MaxkIndex = " << MaxkIndex << endl;
  kVecs.resize(numVecs);
  MagK.resize(numVecs);
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
	  MagK(numVecs) = sqrt(dot(k,k));
	  numVecs++;
	}
      }
    }
  }
}

// void PathClas::Calc_Ylm(int slice,int species)
// {


// }

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
  //  perr<<"Beginning the calcrhok stuff"<<endl;
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
  //  perr<<"ending the calcrhok stuff"<<endl;
}
      
void 
PathClass::UpdateRho_ks(int slice1, int slice2,
			const Array<int,1> &changedParticles, int level)
{
  int skip = 1<<level;
  /// First, reset new to old;
  for (int slice=slice1; slice<=slice2; slice+=skip)
    for (int species=0; species<NumSpecies(); species++)
      for (int ki=0; ki<kIndices.size(); ki++)
	Rho_k[NEWMODE](slice,species,ki) = Rho_k[OLDMODE](slice,species,ki);
  
  /// Now update the new depending on the changed positions.
  for (int pi=0; pi<changedParticles.size(); pi++) {
    int ptcl = changedParticles(pi);
    int species = ParticleSpeciesNum(ptcl);
    for (int slice=slice1; slice <= slice2; slice+=skip) {
      dVec &rnew = Path[NEWMODE](slice,ptcl);
      dVec &rold = Path[OLDMODE](slice,ptcl);
      for (int dim=0;dim<NDIM;dim++){
	complex<double> tempC;
	double phi = rnew[dim] * kBox[dim];
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
	Rho_k[NEWMODE](slice,species,ki) += 
	  C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
	//#ifdef TWO_D
#if NDIM==2
	Rho_k[NEWMODE](slice,species,ki) += C[0](kIndex[0])*C[1](kIndex[1]);
#endif
      }
      for (int dim=0;dim<NDIM;dim++){
	complex<double> tempC;
	double phi = rold[dim] * kBox[dim];
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
	Rho_k[NEWMODE](slice,species,ki) -= 
	  C[0](kIndex[0])*C[1](kIndex[1])*C[2](kIndex[2]);
#endif
	//#ifdef TWO_D
#if NDIM==2
	Rho_k[NEWMODE](slice,species,ki) -= C[0](kIndex[0])*C[1](kIndex[1]);
#endif
      }
    }
  }
}

void 
PathClass::UpdateRho_ks()
{
  ModeType mode = GetMode();
  SetMode(OLDMODE);
  for (int slice=0; slice<NumTimeSlices(); slice++)
    for (int species=0; species<NumSpecies(); species++)
      CalcRho_ks_Fast(slice,species);
  Rho_k[1] = Rho_k[0];
//   SetMode(NEWMODE);
//   for (int slice=0; slice<NumTimeSlices(); slice++)
//     for (int species=0; species<NumSpecies(); species++)
//       CalcRho_ks_Fast(slice,species);
  SetMode(mode);
}



void PathClass::MoveJoin(int oldJoin, int newJoin)
{
  //  perr<<"My time slices is "<<NumTimeSlices()<<endl;
  //  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    //    perr<<Permutation(ptcl)<<endl;
  //  }
  //  perr<<oldJoin<<" "<<newJoin<<" "<<endl;
  //  perr<<"Starting"<<endl;
  //  perr<<Path(OpenLink,OpenPtcl)<<endl;
  bool swappedAlready=false;
  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	Path[OLDMODE](timeSlice,ptcl)=Path[NEWMODE](timeSlice,Permutation(ptcl));
	if (OpenPaths) {
	  if (timeSlice==(int)OpenLink && Permutation(ptcl)==(int)OpenPtcl && !swappedAlready){
	    OpenPtcl[OLDMODE]=ptcl;
	    OpenPtcl[NEWMODE]=ptcl;
	    swappedAlready=true;
	  }
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
	if (OpenPaths)
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
      for (int ptcl=0;ptcl<NumParticles();ptcl++)
	Path[NEWMODE](timeSlice,ptcl)=Path[OLDMODE](timeSlice,ptcl);
    }
  }
  //  perr<<Path(OpenLink,OpenPtcl)<<endl;
  //  perr<<"Ending"<<endl;
}



void PathClass::AcceptCopy(int startSlice,int endSlice, 
			   const Array <int,1> &activeParticles)
{
  ExistsCoupling.AcceptCopy();
  Weight.AcceptCopy();
  NowOpen.AcceptCopy();
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[OLDMODE](Range(startSlice, endSlice), ptcl) = 
      Path[NEWMODE](Range(startSlice, endSlice), ptcl);
    Permutation.AcceptCopy(ptcl);
  }
  //I had this accepting the whole permutation.  Not sure why but have commented it out.  
// <<<<<<< .mine
//   for (int ptcl=0;ptcl<NumParticles();ptcl++)
//     Permutation.AcceptCopy(ptcl);
// =======
  Rho_k[OLDMODE](Range(startSlice,endSlice), Range::all(), Range::all()) =
    Rho_k[NEWMODE](Range(startSlice,endSlice), Range::all(), Range::all());


  if (OpenPaths){
    OpenPtcl.AcceptCopy();
    OpenLink.AcceptCopy();

    for (int counter=0;counter<NumTimeSlices();counter++){
      Path[OLDMODE](counter,NumParticles())=Path[NEWMODE](counter,NumParticles());
    }
    for (int ptcl=0;ptcl<NumParticles();ptcl++){
      Path[OLDMODE](OpenLink[NEWMODE],ptcl)=Path[NEWMODE](OpenLink[NEWMODE],ptcl);
      Path[OLDMODE](0,ptcl)=Path[NEWMODE](0,ptcl);
      Permutation.AcceptCopy(ptcl);
    }
    //    Path[OLDMODE](Range(startSlice,endSlice),NumParticles())=
    //      Path[NEWMODE](Range(startSlice,endSlice),NumParticles());
  }
  //  cerr<<"I am binning from "<<startSlice<<" to "<<endSlice<<endl;
  if (OrderN){
    for (int slice=startSlice;slice<=endSlice;slice++){
      for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
	int ptcl=activeParticles(ptclIndex);
	Cell.ReGrid(slice,ptcl);
	//      Cell.BinParticles(slice);
      }
    }
  }
}

void PathClass::RejectCopy(int startSlice,int endSlice, 
				  const Array <int,1> &activeParticles)
{
  ExistsCoupling.RejectCopy();
  Weight.RejectCopy();
  NowOpen.RejectCopy();
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    Path[NEWMODE](Range(startSlice, endSlice), ptcl) = 
      Path[OLDMODE](Range(startSlice, endSlice), ptcl);
    Permutation.RejectCopy(ptcl);
  }
  //For some reason rejecting the entire permutation?
//   for (int ptcl=0;ptcl<NumParticles();ptcl++)
//     Permutation.RejectCopy(ptcl);

  Rho_k[NEWMODE](Range(startSlice,endSlice), Range::all(), Range::all()) =
    Rho_k[OLDMODE](Range(startSlice,endSlice), Range::all(), Range::all());

  if (OpenPaths){
    OpenPtcl.RejectCopy();
    OpenLink.RejectCopy();

    Path[NEWMODE](Range(startSlice,endSlice),NumParticles())=
      Path[OLDMODE](Range(startSlice,endSlice),NumParticles());
    for (int ptcl=0;ptcl<NumParticles();ptcl++){
      Path[NEWMODE](OpenLink[NEWMODE],ptcl)=Path[OLDMODE](OpenLink[NEWMODE],ptcl);
      Path[NEWMODE](0,ptcl)=Path[OLDMODE](0,ptcl);
      Permutation.RejectCopy(ptcl);
    }

  }

}


void PathClass::ShiftData(int slicesToShift)
{

  ShiftPathData(slicesToShift);
  if (LongRange)
    // ShiftRho_kData(slicesToShift);
    UpdateRho_ks();
  OpenLink.AcceptCopy(); ///the open link has changed and you want to accept it
  RefSlice += slicesToShift;
  while (RefSlice >= TotalNumSlices)
    RefSlice -= TotalNumSlices;
  while (RefSlice < 0)
    RefSlice += TotalNumSlices;
  //  perr<<"My ref slice at particle 0 is "<<Path(RefSlice,0)<<endl;
  if (OpenPaths){

    //    perr<<"Here my open link is "<<OpenLink<<endl;
    int openLinkOld=(int)OpenLink;
    int startSliceOpenLinkProc;
    int endSliceOpenLinkProc;
    SliceRange(SliceOwner(RefSlice),startSliceOpenLinkProc,
	       endSliceOpenLinkProc);
    
    OpenLink=RefSlice-startSliceOpenLinkProc;
    if ((int)OpenLink==0){
      OpenLink=NumTimeSlices()-1;
    }
    //  perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    ////    perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<" "<<slicesToShift<<endl;
    //    if (openLinkOld!=(int)OpenLink)
    //      cerr<<"My links are "<<openLinkOld<<" "<<OpenLink<<" "<<slicesToShift<<endl<<Communicator.MyProc()<<" "<<SliceOwner(RefSlice)<<endl;
    //    assert(openLinkOld==(int)OpenLink);
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
  //  perr<<"Slices to shift are "<<slicesToShift<<endl;
  //  perr<<"I'm in shiftpathdata with numparticles "<<NumParticles()+OpenPaths<<"and Openpaths being "<<OpenPaths<<"and slicesToShift is "<<slicesToShift<<endl;
  //  sleep(10);
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
      perr<<"equal to 0"<<endl;
      OpenLink=NumTimeSlices()-1;
    }

    //    OpenLink=RefSlice;
    //    if ((int)OpenLink==0){
    //      OpenLink=NumTimeSlices()-1;
    //    }
    //    perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    //    if (openLinkOld!=(int)OpenLink)
    //      perr<<"My links are "<<openLinkOld<<" "<<OpenLink<<endl;
    //    assert(openLinkOld==(int)OpenLink);
    OpenLink[OLDMODE]=OpenLink[NEWMODE];
  }
  //  perr<<Path(OpenLink,NumParticles())<<endl;
  //  perr<<"I leave shiftpathdata"<<endl;
  //  sleep(10);

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
  
  //  perr << "procWithRefSlice = " << procWithRefSlice << endl;

  // Do broadcast
  Communicator.Broadcast(procWithRefSlice, buffer);

  // Now, all processors have reference slice.  Note that the 
  // labeling of particles may not be right since with haven't
  // propogated the permuation matrices.
  
  // Now, copy into path
  for (int ptcl=0; ptcl<NumParticles(); ptcl++)
    RefPath(ptcl) = buffer(ptcl);
}


// Combines the permutation vectors of all the processors.  Note:
// only processor 0 gets the result.  All others get junk.
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
  if (Communicator.MyProc() == 0) {
    for (int pi=0; pi < numPtcls; pi++) {
      int ptcl = pi;
      for (int proc=0; proc<numProcs; proc++)
	ptcl = permMat(proc, ptcl);
      permVec(pi) = ptcl;
    }
  }
}

///Must start diagonal or this is goign to break because of where
///you are starting the openlink at 0
void PathClass::InitOpenPaths()
{
  perr<<"Starting to initialize"<<endl;
  if (OpenPaths){

    perr<<"openpaths"<<endl;
    SetMode(OLDMODE);

    if (Communicator.MyProc()==0)
      OpenLink=NumTimeSlices()-1;
    else
      OpenLink=-1;
    perr<<"Set up the open link"<<endl;
    OpenPtcl=Species(OpenSpeciesNum).FirstPtcl;
    perr<<"set up the open particle"<<endl;
    SetMode(NEWMODE);

    if (Communicator.MyProc()==0)
      OpenLink=NumTimeSlices()-1;
    else 
      OpenLink=-1;
    OpenPtcl=0;
    if (Communicator.MyProc()==0){
      perr<<"Preparing for moving things in path around"<<endl;
      Path[OLDMODE]((int)OpenLink,NumParticles())=Path[OLDMODE]((int)OpenLink,(int)OpenPtcl);
      perr<<"Moved the first thing"<<endl;
      Path[NEWMODE]((int)OpenLink,NumParticles())=Path[NEWMODE]((int)OpenLink,(int)OpenPtcl);
      perr<<"Moved the second thing"<<endl;
      Path[OLDMODE](0,NumParticles())=Path[OLDMODE](0,(int)OpenPtcl);
      perr<<"Moved the first thing"<<endl;
      Path[NEWMODE](0,NumParticles())=Path[NEWMODE](0,(int)OpenPtcl);
    }
  }
  perr<<"Initialized the open paths"<<endl;
}

void
PathClass::WarpAtoB(dVec &pos)
{
  double weightSum = 0.0;
  for (int i=0; i<IonConfigs[0].size(); i++) {
    dVec disp = pos - IonConfigs[0](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    weightSum += 1.0/(dist2*dist2);
  }
  dVec shift;
  for(int i=0; i<NDIM; i++)
    shift[i] = 0;
  for (int i=0; i<IonConfigs[0].size(); i++) {
    dVec disp = pos - IonConfigs[0](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    double weight = 1.0/(dist2*dist2*weightSum);
    shift = shift + weight*(IonConfigs[1](i)-IonConfigs[0](i));
  }
  pos = pos + shift;
}

void
PathClass::WarpBtoA(dVec &pos)
{
  double weightSum = 0.0;
  for (int i=0; i<IonConfigs[0].size(); i++) {
    dVec disp = pos - IonConfigs[1](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    weightSum += 1.0/(dist2*dist2);
  }
  dVec shift;
  for(int i=0; i<NDIM; i++)
    shift[i] = 0;
  for (int i=0; i<IonConfigs[1].size(); i++) {
    dVec disp = pos - IonConfigs[1](i);
    PutInBox(disp);
    double dist2 = dot(disp, disp);
    double weight = 1.0/(dist2*dist2*weightSum);
    shift = shift + weight*(IonConfigs[0](i)-IonConfigs[1](i));
  }
  pos = pos + shift;
}

void 
PathClass::WarpPaths (int ionSpecies)
{
  SpeciesClass &ions = Species(ionSpecies);
  int N = ions.NumParticles;
  Array<dVec,1> Rold(N), Rnew(N), Rdelta(N);
  for (int i=0; i<N; i++) {
    Rold(i) = Path[0](0,i+ions.FirstPtcl);
    Rnew(i) = Path[1](0,i+ions.FirstPtcl);
    Rdelta(i) = Rnew(i)-Rold(i);
    PutInBox(Rdelta(i));
  }
  SetMode(OLDMODE);
  for (int elec=0; elec<NumParticles(); elec++) {
    if (ParticleSpeciesNum(elec) != ionSpecies) {
      for (int slice=0; slice<NumTimeSlices(); slice++) {
	dVec disp;
	double dist;
	double totalWeight = 0.0;
	for (int ion=0; ion<N; ion++) {
	  DistDisp(slice, elec, ion+ions.FirstPtcl, dist, disp);
	  totalWeight += 1.0/(dist*dist*dist*dist);
	}
	for (int ion=0; ion<N; ion++) {
	  DistDisp(slice, elec, ion+ions.FirstPtcl, dist, disp);
	  double weight = 1.0/(dist*dist*dist*dist*totalWeight);
	  Path[1](slice, elec) += weight*Rdelta(ion);
	}
	Path[0](slice, elec) = Path[1](slice, elec);
      }
    }
  }
}
