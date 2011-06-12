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

#include "VacancyDensity.h"


int VacancyDensityClass::IntoGrid(double num, int dim)
{
  dVec box=PathData.Path.GetBox();
  //  double maxLen=max(box[0],max(box[1],box[2]));
  double lent=box[dim];
  while (num>lent/2) 
    num-=lent;
  while (num<-lent/2)
    num+=lent;
  num+=lent/2.0;
  int myNum=(int)floor(num/(lent/(double)(Grid.extent(dim))));
  return myNum;
}


// ///Only works in cubic box
// int VacancyDensityClass::IntoGrid(double num)
// {
//   dVec box=PathData.Path.GetBox();
//   double maxLen=max(box[0],max(box[1],box[2]));
//   while (num>maxLen/2) 
//     num-=maxLen;
//   while (num<-maxLen/2)
//     num+=maxLen;
//   num+=maxLen/2.0;
//   int myNum=(int)floor(num/(maxLen/(double)(Grid.extent(0))));
//   //  //  cerr<<num;
//   //  cerr<<maxLen/(double)(Grid.extent(0));
//   //  cerr<<"My num is "<<myNum<<endl;
//   return myNum;
// }


// Fix to include final link between link M and 0
void VacancyDensityClass::Accumulate()
{
  for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++){
    VacancyHelper.FindVacancy(slice);
    //    cerr<<"The vacancy is located at "<<VacancyHelper.VacancyArray(0)<<" "<<slice<<endl;
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      NumSamples++;

      dVec pos=VacancyHelper.FixedLoc(VacancyHelper.VacancyArray(0))-
	PathData.Path(slice,ptcl);
      int nx=IntoGrid(pos[0],0);
      int ny=IntoGrid(pos[1],1);
      int nz=IntoGrid(pos[2],2);
      assert(nx<Grid.extent(0));
      assert(ny<Grid.extent(1));
      assert(nz<Grid.extent(2));
      Grid(nx,ny,nz)=Grid(nx,ny,nz)+1;
    }
	
  }
}

void VacancyDensityClass::WriteBlock()
{
  double norm = 1.0/((double)NumSamples);
  Grid=Grid/norm;
  GridVar.Write(Grid);
  GridVar.Flush();
  NumSamples = 0;
  Grid=0.0;

}

void VacancyDensityClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);

  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Grid");
  }
  VacancyHelper.Read(in);
  
}



void VacancyDensityClass::WriteInfo()
{


}
