#include <sstream>

#include <Common/IO/IO.h>
#include <blitz/array.h>
#include "Observables/ObservableBase.h"
 #define FORT(name) name ## _
 #define F77_LSAPR F77_FUNC(lsapr,LSAPR)

 using namespace blitz;

string IntToString ( int number )
{
  ostringstream oss;

  // Works just like cout
  oss<< "ptcl."<<number;

  // Return the underlying string
  return oss.str();
}

void 
PutInBox (dVec &v,Array<double,1> Box)
{
  for (int i=0; i<NDIM; i++) {
    double n = -floor(v(i)*(1.0/Box(i))+0.5);
    v(i) += n*Box(i);
  }
}

int closestLatticeSite(Array<double,2> positions,double xpos, double ypos, double zpos,Array<double,1> box)
{
  int closestPosition=0;
  double closestDist=99999;
  for (int i=0;i<positions.extent(0);i++){
    dVec diff;
    diff[0]=positions(i,0)-xpos;
    diff[1]=positions(i,1)-ypos;
    diff[2]=positions(i,2)-zpos;
    PutInBox(diff,box);
    double dist=dot(diff,diff);
    if (dist<closestDist){
      closestDist=dist;
      closestPosition=i;
    }
  }
  return closestPosition;
}

double EvaluatePerm(Array<int,1> Perm,
		  Array<double,2> positions,
		  Array<double,4> oldPaths,
		  int mcStep,
		  int slice,
		  Array<double,1> box)
{
  double totalVal=0.0;
  for (int i=1;i<Perm.size();i++){
    dVec diff;
    diff[0]=oldPaths(mcStep,i-1,slice,0)-positions(Perm(i)-1,0);
    diff[1]=oldPaths(mcStep,i-1,slice,1)-positions(Perm(i)-1,1);
    diff[2]=oldPaths(mcStep,i-1,slice,2)-positions(Perm(i)-1,2);
    PutInBox(diff,box);
    double dist2=dot(diff,diff);
    totalVal+=max(dist2,0.25*3.67*3.67);
  }
  return totalVal;


}
extern "C" void 
F77_LSAPR (int *n,double* c, int *perm);

int main()
{

  Array<double,3> pathArray;
  Array<double,2> positions;
  
  IOSectionClass InputFile;
  InputFile.OpenFile("pimc.in");
  InputFile.OpenSection("Observables");
  string outFileBase;
  InputFile.ReadVar("OutFileBase",outFileBase);
  cerr<<outFileBase<<endl;
  int num=0;
  InputFile.OpenSection("Observable",num);
  string typeString;
  InputFile.ReadVar("type",typeString);
  cerr<<typeString<<endl;
  bool hasVariable=InputFile.ReadVar("LocationsToCompare",positions);
  InputFile.CloseSection();
  while (!hasVariable){
    num++;
    InputFile.OpenSection("Observable",num);
    InputFile.ReadVar("type",typeString);
    cerr<<typeString<<endl;
    hasVariable=InputFile.ReadVar("LocationsToCompare",positions);
    InputFile.CloseSection();
  }
  InputFile.CloseSection();
  Array <double,1> box;
  InputFile.OpenSection("System");
  InputFile.ReadVar("Box",box);
  InputFile.CloseSection();
  InputFile.CloseFile();
  cerr<<"My box is "<<box<<endl;
  ofstream positionsFile;
  positionsFile.open("positions");
  
  for (int i=0;i<positions.size();i++){
    dVec pos;
    pos(0)=positions(i,0); pos(1)=positions(i,1); pos(2)=positions(i,2);
    PutInBox(pos,box);
    positionsFile << pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
  }

  IOSectionClass IOS; 
  IOS.OpenFile("SingleVacancyLowT2.2.h5");
  IOS.OpenSection("Observables");
  IOS.OpenSection("PathDump");
  Array<double,4> oldPaths;
  assert(IOS.ReadVar("Path",oldPaths));
  IOS.CloseSection();
  IOS.CloseSection();
  IOS.CloseFile();

  cerr << "My paths are of size"  << oldPaths.extent(0) << " "
	 << oldPaths.extent(1)<<" " << oldPaths.extent(2) << " "
	 << oldPaths.extent(3)<<      endl;
    


  bool firstTime=true;
  int numMCTime=oldPaths.extent(0);
  int numSlice=oldPaths.extent(2);
  int numPtcl=oldPaths.extent(1);
  int numDim=oldPaths.extent(3);

  IOSectionClass outFile;
  outFile.NewFile("toVisualize.h5");
  outFile.NewSection("System");
  outFile.WriteVar("Box",box);
  outFile.WriteVar("NumTimeSlices",numSlice);
  outFile.NewSection("Species");
  outFile.WriteVar("lambda",6.05);
  outFile.WriteVar("Name","Vacancy");
  outFile.WriteVar("NumParticles",1);
  outFile.CloseSection();
  outFile.CloseSection();
  outFile.NewSection("Observables");
  outFile.NewSection("PathDump");

  Array<int,1> permVec;
  permVec.resize(1);
  for (int counter=0;counter<permVec.size();counter++)
    permVec(counter)=counter;

  int numEmptySites=1;
  pathArray.resize(1,numSlice,numDim);
  Array<double,2> DistTable(1,1,ColumnMajorArray<2>());
  DistTable.resize(numPtcl+numEmptySites, numPtcl+numEmptySites);
  Array<int,1> Perm;
  Perm.resize(numPtcl+numEmptySites);
  DistTable=0.0;
  cerr<<"Positions extent 0 is "<<positions.extent(0)<<endl;
  IOVarBase *PathVar;
  IOVarBase *PermVar;
  int prePerm=0;
  int nextPerm=0;
  ////HACK! BUG!  for (int mcStep=0;mcStep<numMCTime;mcStep++){
  ofstream outfiles[numPtcl+numEmptySites];
  for (int i=0;i<numPtcl+numEmptySites;i++){
    string fileName=IntToString(i);
    outfiles[i].open(fileName.c_str());
  }
  for (int mcStep=numMCTime-1;mcStep<numMCTime;mcStep++){
  for (int slice=0;slice<numSlice;slice++){
//     if (slice>=1){
//       cerr<<"My perm value is "<<slice<<" "<<EvaluatePerm(Perm,
// 							  positions,
// 							  oldPaths,
// 							  mcStep,
// 							  slice,
// 							  box)<<endl;

//     }
    //    DistTable=0.0;
    for (int latticeSite=0;latticeSite<positions.extent(0);latticeSite++){
      for (int ptcl=0;ptcl<numPtcl;ptcl++){
	dVec disp;
	double dist2;
	disp(0)=oldPaths(mcStep,ptcl,slice,0)-positions(latticeSite,0);
	disp(1)=oldPaths(mcStep,ptcl,slice,1)-positions(latticeSite,1);
	disp(2)=oldPaths(mcStep,ptcl,slice,2)-positions(latticeSite,2);
	PutInBox(disp,box);
	dist2=dot(disp,disp);
	//	DistTable(latticeSite,ptcl+numEmptySites)=max(dist2,0.25*3.67*3.67);
	DistTable(latticeSite,ptcl+numEmptySites)=dist2;
      }
    }
    int n =positions.extent(0);
    //    Perm=0;
    F77_LSAPR (&n,DistTable.data(),Perm.data());
    prePerm=nextPerm;
    nextPerm=Perm(0);
    cerr<<"My CORRECT perm value is "<<slice<<" "<<Perm(0)<<" "
	<<EvaluatePerm(Perm,
		       positions,
		       oldPaths,
		       mcStep,
		       slice,
		       box)<<endl;
    cerr<<"CORRECT: "<<Perm<<endl;
    for (int ptcl=0;ptcl<Perm.size();ptcl++)
      outfiles[ptcl] << positions(Perm(ptcl)-1,0)<<" "<<positions(Perm(ptcl)-1,1)<<" "<<positions(Perm(ptcl)-1,2)<<endl;
    
//       ///HACK!      

//     //    DistTable=0.0;
//     for (int latticeSite=0;latticeSite<positions.extent(0);latticeSite++){
//       for (int ptcl=0;ptcl<numPtcl;ptcl++){
// 	dVec disp;
// 	double dist2;
// 	disp(0)=oldPaths(mcStep,ptcl,slice,0)-positions(latticeSite,0);
// 	disp(1)=oldPaths(mcStep,ptcl,slice,1)-positions(latticeSite,1);
// 	disp(2)=oldPaths(mcStep,ptcl,slice,2)-positions(latticeSite,2);
// 	PutInBox(disp,box);
// 	dist2=dot(disp,disp);
// 	DistTable(latticeSite,ptcl+numEmptySites)=dist2;
//       }
//     }

//     cerr<<"Perm A"<<Perm<<endl;

// //     for (int k=1;k<DistTable.extent(1);k++)
// //       DistTable(prePerm-1,k)=5*20*20*20.0;
//       //      DistTable(k,Perm(0)-1)=5*20*20*20.0;
//     //      Perm=0;
//       n =positions.extent(0);

//       F77_LSAPR (&n,DistTable.data(),Perm.data());
//       cerr<<"My perm value afteris "<<slice<<" "<<Perm(0)<<" "
// 	  <<EvaluatePerm(Perm,
// 			 positions,
// 			 oldPaths,
// 			 mcStep,
// 			 slice,
// 			 box)<<endl;

//       cerr<<"Perm B"<<Perm<<endl;

//       //HACK!

    pathArray(0,slice,0)=positions(Perm(0)-1,0);
    pathArray(0,slice,1)=positions(Perm(0)-1,1);
    pathArray(0,slice,2)=positions(Perm(0)-1,2);
    dVec distCheck;
    distCheck[0]=pathArray(0,(slice+numSlice-1) % numSlice,0)-
      pathArray(0,slice,0);
    distCheck[1]=pathArray(0,(slice+numSlice-1) % numSlice,1)-
      pathArray(0,slice,1);
    distCheck[2]=pathArray(0,(slice+numSlice-1) % numSlice,2)-
      pathArray(0,slice,2);
    PutInBox(distCheck,box);
    double dist=sqrt(dot(distCheck,distCheck));
    cerr<<"Data: "<<mcStep<<" "
	<<pathArray(0,slice,0)<<" "
	<<pathArray(0,slice,1)<<" "
	<<pathArray(0,slice,2)<<" "
	<<dist<<endl;
    for (int i=0;i<Perm.size();i++)
      cerr<<Perm(i)<<" ";
    cerr<<endl;
    cerr<<"Starting closest"<<endl;
    for (int i=0;i<oldPaths.extent(1);i++)
      cerr<<closestLatticeSite(positions,oldPaths(mcStep,i,slice,0),oldPaths(mcStep,i,slice,1),oldPaths(mcStep,i,slice,2),box)+1<<" ";
    cerr<<endl;
    cerr<<"Ending closest"<<endl;
    cerr<<endl;
    cerr<<oldPaths(mcStep,0,slice,0)<<" "
	<<oldPaths(mcStep,0,slice,1)<<" "
	<<oldPaths(mcStep,0,slice,2)<<" "<<endl;
    //    cerr<<oldPaths(mcStep,slice,150,0)<<" "
    //	<<oldPaths(mcStep,slice,150,1)<<" "
    //	<<oldPaths(mcStep,slice,150,2)<<" "<<endl;

    //    cerr<<"Perm: "<<Perm<<endl;

      //    cerr<<Perm(0)<<endl;
  }
  cerr<<"About to write"<<endl;
  cerr<<pathArray<<endl;
   if (firstTime){
     cerr<<"the perm vec size is "<<permVec.extent(0)<<endl;
     Array<int,2> tensor2;
     tensor2.resize(1,permVec.extent(0));
     tensor2(0,Range::all())=permVec;
     outFile.WriteVar("Permutation",tensor2);
     PermVar = outFile.GetVarPtr("Permutation");

     Array<double,4> tensor;
     tensor.resize(1,pathArray.extent(0), pathArray.extent(1), pathArray.extent(2));
     tensor(0,Range::all(),Range::all(),Range::all()) = pathArray;
     outFile.WriteVar ("Path", tensor);
     PathVar = outFile.GetVarPtr("Path");
     firstTime=false;
   }
   else{
     PathVar->Append(pathArray);
     PermVar->Append(permVec);
   }
  }
   cerr<<"Written"<<endl;
   outFile.CloseSection();
//   cerr<< "A"<<endl;
   outFile.CloseSection();
   outFile.CloseFile();
//   cerr<<"B"<<endl;
//   //  //  outFile.CloseFile();
//   //  cerr<<"C"<<endl;
  
  
}
