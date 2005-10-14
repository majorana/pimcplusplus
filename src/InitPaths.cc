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
