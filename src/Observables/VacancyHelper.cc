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


#include "VacancyHelper.h"
#include <Common/Blitz.h>
#include <blitz/array.h>
#include <Common/IO/IO.h>
#include "../Common.h"

#define FORT(name) name ## _
#define F77_LSAPR F77_FUNC(lsapr,LSAPR)

using namespace blitz;



extern "C" void 
F77_LSAPR (int *n,double* c, int *perm);

void 
VacancyHelperClass::FindVacancy(int slice)
{
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  Array<double,2> DistTable(1,1,ColumnMajorArray<2>());
  DistTable.resize(FixedLoc.size(), FixedLoc.size());
  Array<int,1> Perm;
  Perm.resize(FixedLoc.size());
  DistTable=0.0;
  int numEmptySites=FixedLoc.size()-PathData.Path.NumParticles();
  for (int latticeSite=0;latticeSite<FixedLoc.size();latticeSite++){
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      dVec disp;
      double dist2;
      disp=PathData.Path(slice,ptcl)-FixedLoc(latticeSite);
      PathData.Path.PutInBox(disp);
      dist2=dot(disp,disp);
      DistTable(latticeSite,ptcl+numEmptySites)=dist2;
    }
    int n =FixedLoc.size();
    F77_LSAPR (&n,DistTable.data(),Perm.data());
    for (int i=0;i<NumVacancies;i++)
      //subtract by 1 is to compensate for the fact that fortran is
      //indexing by 1 
      VacancyArray(i)=Perm(i)-1;
  }
}    


void VacancyHelperClass::Read(IOSectionClass &in)
{  

  Array<double,2> positions;
  assert(in.ReadVar("LatticeSites",positions));
  FixedLoc.resize(positions.extent(0));
  ///Verify you used the right number of points to compare against
  assert(positions.extent(1)==NDIM);
  dVec pos;
  for (int loc=0;loc<FixedLoc.size(); loc++){
    for (int dim=0; dim<NDIM; dim++)
      pos(dim) = positions(loc,dim);
    FixedLoc(loc) = pos;
  }      
  NumVacancies=FixedLoc.size()-PathData.Path.NumParticles();
  cerr<<"Size of Vacancy ARray is "<<NumVacancies<<endl;
  VacancyArray.resize(NumVacancies);
  
}




