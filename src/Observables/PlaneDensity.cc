#include "PlaneDensity.h"

///Only works in cubic box
int PlaneDensityClass::IntoGrid(double num)
{
  dVec box=PathData.Path.GetBox();
  double maxLen=max(box[0],max(box[1],box[2]));
  while (num>maxLen/2) 
    num-=maxLen;
  while (num<-maxLen/2)
    num+=maxLen;
  num+=maxLen/2.0;
  int myNum=(int)floor(num/(maxLen/(double)(Grid.extent(0))));
  //  //  cerr<<num;
  //  cerr<<maxLen/(double)(Grid.extent(0));
  //  cerr<<"My num is "<<myNum<<endl;
  return myNum;
}


// Fix to include final link between link M and 0
void PlaneDensityClass::Accumulate()
{
  for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
      if (abs(PathData.Path(slice,ptcl)[2])<1.85){
	NumSamples++;
	int nx=IntoGrid(PathData.Path(slice,ptcl)[0]);
	int ny=IntoGrid(PathData.Path(slice,ptcl)[1]);
	if (nx<NumSamples && ny<NumSamples)
	  Grid(nx,ny)=Grid(nx,ny)+1;
      }

}

void PlaneDensityClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  Grid=Grid/norm;
  GridVar.Write(Grid);
  GridVar.Flush();
  NumSamples = 0;
  Grid=0.0;

}

void PlaneDensityClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);

  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Grid");
  }
  
  
}



void PlaneDensityClass::WriteInfo()
{


}
