#include "PathClass.h"
#include "Actions/ActionsClass.h"


void PathClass::ReadSqueeze(string fileName,bool replicate)
{
// This was modified on Jan 19 2005 to read in a set of classical (P=1) configs and duplicate them to produce a set of PIMC (P>1) configs.  -jg
  IOSectionClass inFile;
  assert (inFile.OpenFile(fileName.c_str()));
  inFile.OpenSection("System");
  Array<double,3> oldBox;
  inFile.ReadVar("Box",oldBox);
  inFile.CloseSection();
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<double,4> oldPaths; //(58,2560,2,3);
  assert(inFile.ReadVar("Permutations",Permutation.data()));
  assert(inFile.ReadVar("Path",oldPaths));
  perr << "My paths are of size"  << oldPaths.extent(0) << " "
       << oldPaths.extent(1)<<" " << oldPaths.extent(2) << endl;
  
  for (int ptcl=0;ptcl<NumParticles();ptcl++){
    for (int slice=0; slice<NumTimeSlices(); slice++) {
      dVec pos;
      pos = 0.0;
      for (int dim=0; dim<NDIM; dim++)
	if (replicate){
          pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,0,dim)*(Box[dim]/oldBox(dim));
        }
        else{
	  pos(dim) = oldPaths(oldPaths.extent(0)-1,ptcl,slice,dim)*(Box[dim]/oldBox(dim));
        }
      cerr<<"I'm putting the slice "<<slice<<" and the ptcl "<<ptcl<<"as "<<Path(slice,ptcl)<<endl;
      Path(slice,ptcl) = pos;
    }      
  }
  inFile.CloseSection();
  inFile.CloseSection();
  inFile.CloseFile();
}

void PathClass::ReadOld(string fileName,bool replicate)
{
// This was modified on Jan 19 2005 to read in a set of classical (P=1) configs and duplicate them to produce a set of PIMC (P>1) configs.  -jg
  cerr<<"Trying to read old"<<endl;
  IOSectionClass inFile;
  assert (inFile.OpenFile(fileName.c_str()));
  inFile.OpenSection("Observables");
  inFile.OpenSection("PathDump");
  Array<double,4> oldPaths; //(58,2560,2,3);
  
  assert(inFile.ReadVar("Path",oldPaths));
  perr << "My paths are of size"  << oldPaths.extent(0) << " "
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
  //  cerr<<"Hello"<<endl;
  //  assert(1==2);
  assert(in.OpenSection ("Particles"));
  int numSpecies = NumSpecies();
  // Now initialize the Path
  for (int speciesIndex=0; speciesIndex<numSpecies; speciesIndex++)
  {
    SpeciesClass &species = *SpeciesArray(speciesIndex);
    assert(in.OpenSection("Species", speciesIndex));
    string InitPaths;
    cerr<<"about to read the string"<<endl;
    in.ReadVar ("InitPaths", InitPaths);
    cerr<<"Read "<<InitPaths<<endl;
    string Replicate;
    in.ReadVar ("Replicate", Replicate);
    if (InitPaths == "RANDOM") {
      perr << "Don't know how to do RANDOM yet.\n";
      exit(1);
    }
    else if (InitPaths == "CUBIC") {
      int num = species.NumParticles;
      bool isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
      if (!isCubic) {
	perr << "A cubic box is current required for cubic initilization\n";
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
    else if(InitPaths=="UniformSphere") {
      dVec r0,r;
      double SphereRadius;
      assert(in.ReadVar("SphereRadius",SphereRadius));
      //      double SphereRadius=31;
      SpeciesClass &species = *SpeciesArray(0);
      double sigma=sqrt(2*species.lambda*tau);
      for (int ptcl =0;ptcl<NumParticles();ptcl++){
        Random.LocalGaussianVec(1,r0);//under common
        cerr<<"in PathClass"<<r0*SphereRadius/sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2])<<endl;
	for (int slice=0;slice <NumTimeSlices();slice++){
          SetPos(slice,ptcl,r0*SphereRadius/sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]));                      
	}
      }
    }
    else if(InitPaths=="UniformPointsSphere") {
      dVec r0,r;
      SpeciesClass &species = *SpeciesArray(0);
      double SphereRadius;
      assert(in.ReadVar("SphereRadius",SphereRadius));
      //      double SphereRadius=31;
      double sigma=sqrt(2*species.lambda*tau);
      for (int ptcl =0;ptcl<NumParticles();ptcl++){
        Random.LocalGaussianVec(1,r0);//under common
        for (int slice=0;slice <NumTimeSlices();slice++){
          Random.LocalGaussianVec(sigma,r);
          r=r0+r;
          SetPos(slice,ptcl,r*SphereRadius/sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));
                                                                                
        }
      }
    }
    else if (InitPaths == "BCC") {
      int num = species.NumParticles;
      bool isCubic = (Box[0]==Box[1]) && (Box[1]==Box[2]);
      if (!isCubic) {
	perr << "A cubic box is current required for cubic initilization\n";
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
// 	fprintf (stderr, "BCC ptcl %d position = [%8.4f %8.4f %8.4f]\n",
// 		 ptcl, r[0], r[1], r[2]);
	for (int slice=0; slice<NumTimeSlices(); slice++) 
	  Path(slice,ptcl) = r;
      }
    }
    else if (InitPaths=="ALLFIXED"){
      Array<double,3> Positions;
      assert(in.ReadVar("Positions",Positions));
      cerr<<"My time slices are "<<NumTimeSlices()<<endl;
      assert(Positions.extent(0)==NumTimeSlices()-1);
      assert(Positions.extent(1)==species.NumParticles);
      for (int ptcl=species.FirstPtcl;
	   ptcl<=species.LastPtcl;ptcl++){
	for (int slice=0;slice<NumTimeSlices()-1;slice++){
	  dVec pos;
	  pos=0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(slice,ptcl-species.FirstPtcl,dim);
	  Path(slice,ptcl) = pos;	  
	}
	dVec pos;
	for (int dim=0; dim<species.NumDim; dim++)
	  pos(dim) = Positions(0,ptcl-species.FirstPtcl,dim);
	Path(NumTimeSlices()-1,ptcl)=pos;
      }
    }
    else if (InitPaths == "FIXED") {
      Array<double,2> Positions;
      assert (in.ReadVar ("Positions", Positions));
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      for (int ptcl=species.FirstPtcl; 
	   ptcl<=species.LastPtcl; ptcl++){
	for (int slice=0; slice<=NumTimeSlices(); slice++) {
	  dVec pos;
	  pos = 0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(ptcl-species.FirstPtcl,dim);
	  Path(slice,ptcl) = pos;
	}      
      }
    }    
    else if (InitPaths == "ALLPATHS") {
      Array<double,3> Positions;
      assert (in.ReadVar ("Positions", Positions));
      perr<<"NumTimeSlices: "<<TotalNumSlices<<" "<<Positions.extent(0);
      assert (Positions.extent(0) == species.NumParticles);
      assert (Positions.extent(1) == TotalNumSlices);
      assert (Positions.extent(2) == species.NumDim);
      for (int ptcl=species.FirstPtcl; 
	   ptcl<=species.LastPtcl; ptcl++){
	for (int slice=0; slice<TotalNumSlices; slice++) {
	  perr<<ptcl;
	  dVec pos;
	  pos = 0.0;
	  for (int dim=0; dim<species.NumDim; dim++)
	    pos(dim) = Positions(ptcl-species.FirstPtcl,slice,dim);
	  Path(slice,ptcl) = pos;
	}      
	int slice=NumTimeSlices()-1;
	dVec pos;
	pos = 0.0;
	for (int dim=0; dim<species.NumDim; dim++)
	  pos(dim) = Positions(ptcl-species.FirstPtcl,0,dim);
	Path(slice,ptcl) = pos;
	
	
      }
    }    
    else if (InitPaths == "FILE"){
      cerr<<"I'm going to read the file now"<<endl;
      bool replicate = false;
      if (Replicate == "ON"){
        replicate = true;
        perr << "Replicate 'ON'; Using time slice 0 for all time slices." << endl;
      }
      string pathFile;
      assert(in.ReadVar("File",pathFile));
      ReadSqueeze(pathFile,replicate);
    }
    else if (InitPaths == "SQUEEZE"){
      string pathFile;
      assert(in.ReadVar("File",pathFile));
    }
    else if (InitPaths == "LEVIFLIGHT") {
      Array<double,2> Positions;
      assert (in.ReadVar ("Positions", Positions));
      assert (Positions.rows() == species.NumParticles);
      assert (Positions.cols() == species.NumDim);
      Array<dVec,1> R0(species.NumParticles);
      for (int ptcl=0; ptcl<species.NumParticles; ptcl++) 
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
      perr << "Unrecognize initialization strategy " 
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
//  perr << "PRINT MolRef corrected" << endl;
  for(int m = 0;m < MolRef.size(); m++){
    int ref = MolRef(m);
    while(ref >= numMol){
      ref -= numMol;
    MolRef(m) = ref;
    }
    perr << m << " " << MolRef(m) << endl;
  } 
  string openSpeciesName;
  if (OpenPaths) {
    assert(in.ReadVar("OpenSpecies",openSpeciesName));
    OpenSpeciesNum=SpeciesNum(openSpeciesName);
    InitOpenPaths();
  }  
  if (OrderN){
    Array<int,1>  numGrid(3);
    //  numGrid(0)=6*2;
    //  numGrid(1)=6*2;
    //  numGrid(2)=5*2;
    numGrid(0)=10;
    numGrid(1)=10;
    numGrid(2)=10;
    
    //  Cell=new GridClass(*this);
    Cell.Init(Box,numGrid);
    //  Cell.BuildNeighborGrids();
    //  cerr<<"I am now printing neighbor grids"<<endl;
    //  Cell.PrintNeighborGrids();
    for (int slice=0;slice<NumTimeSlices();slice++)
      Cell.BinParticles(slice);
    //  cerr<<"I have binned them"<<endl;
    //Cell.PrintParticles(0);
  }
  //Everything needs to be accepted
  Path.AcceptCopy();
  Permutation.AcceptCopy();
  Rho_k.AcceptCopy();
  BroadcastRefPath();
  RefPath.AcceptCopy();
  Weight.AcceptCopy();
  ExistsCoupling.AcceptCopy();
}
