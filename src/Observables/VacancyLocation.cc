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


#include "VacancyLocation.h"
#include <blitz/array.h>

#define FORT(name) name ## _
#define F77_LSAPR F77_FUNC(lsapr,LSAPR)

using namespace blitz;



extern "C" void 
F77_LSAPR (int *n,double* c, int *perm);



// void VacancyLocClass::PrintNearbySites()
// {

//   for (int site=0;site<FixedLoc.size();site++){
//     list<int>::iterator neighborIter;
//     cerr<<"Site: "<<site<<endl;
//     for (neighborIter=Neighbors[site].begin();
// 	 neighborIter!=Neighbors[site].end();
// 	 neighborIter++){
//       cerr<<(*neighborIter)<<", ";
//     }
//     cerr<<endl;
//   }
// }


// void VacancyLocClass::TabulateNearbySites()
// {
//   double minDist=10*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
//   for (int siteA=0;siteA<FixedLoc.size();siteA++){
//     for (int siteB=0;siteB<FixedLoc.size();siteB++){
//       double currDist=
// 	PathData.Path.MinImageDistance(FixedLoc(siteA),FixedLoc(siteB));
//       if (currDist<minDist && abs(currDist)>0.01)
// 	minDist=currDist;
//     }
//   }
//   double epsilon=0.1;
//   for (int siteA=0;siteA<FixedLoc.size();siteA++){
//     for (int siteB=0;siteB<FixedLoc.size();siteB++){
//       double currDist=
// 	PathData.Path.MinImageDistance(FixedLoc(siteA),FixedLoc(siteB));
//       if (currDist<minDist+epsilon && abs(currDist)>0.01)
// 	Neighbors[siteA].push_back(siteB);
//     }
//   }
// }

// bool VacancyLocClass::NeighborsVacancyFree(int site)
// {
//   list<int>::iterator neighborIter;
//   for (neighborIter=Neighbors[site].begin();
//        neighborIter!=Neighbors[site].end();
//        neighborIter++){
//     if (TempLoc(*neighborIter)==0)
//       return false;
//   }
//   return true;
// }

// bool VacancyLocClass::NeighborsNotDoublyOccupied(int site)
// {
//   list<int>::iterator neighborIter;
//   for (neighborIter=Neighbors[site].begin();
//        neighborIter!=Neighbors[site].end();
//        neighborIter++){
//     if (TempLoc(*neighborIter)>1)
//       return false;
//   }
//   return true;
// }

// void VacancyLocClass::ReturnDispFromASite(dVec loc1,dVec loc2)
// {
//   for (int loc=0;loc<FixedLoc.size();loc++){
//     dVec dispCompare=FixedLoc(loc)-FixedLoc(0);
//     PathData.Path.PutInBox(dispCompare);
//     if (dot(disp1-dispCompare,disp1-dispCompare)<=1e-5){
//       disp1Found=true;
//       HistogramDisp(loc)=HistogramDisp(loc)+1.0;
//     }
//     if (dot(disp2-dispCompare,disp2-dispCompare)<=1e-5){
//       HistogramDisp(loc)=HistogramDisp(loc)+1.0;
//       disp2Found=true;
//     }
    
// 	}
// }

void VacancyLocClass::Accumulate()
{
  TimesCalled++;
//   if (TimesCalled % DumpFreq==0)
//     WriteBlock();

//   if ((TimesCalled % Freq)!=0){
//     return;
//   }

  if ((TimesCalled % Freq)!=0){
    return;
  }
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
 
    Array<double,2> DistTable(1,1,ColumnMajorArray<2>());
    DistTable.resize(FixedLoc.size(), FixedLoc.size());
    Array<int,1> Perm;
    Perm.resize(FixedLoc.size());
    DistTable=0.0;
    int numEmptySites=FixedLoc.size()-PathData.Path.NumParticles();

    for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++){
      for (int latticeSite=0;latticeSite<FixedLoc.size();latticeSite++){
	for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
	dVec disp;
	double dist2;
	disp=PathData.Path(slice,ptcl)-FixedLoc(latticeSite);
	PathData.Path.PutInBox(disp);
	dist2=dot(disp,disp);
	DistTable(latticeSite,ptcl+numEmptySites)=dist2;
	}
      }
      int n =FixedLoc.size();
      F77_LSAPR (&n,DistTable.data(),Perm.data());
      if (numEmptySites==2){
	int vacancy1=Perm(0)-1;
	int vacancy2=Perm(1)-1;
	dVec disp1=FixedLoc(vacancy1)-FixedLoc(vacancy2);
	dVec disp2=FixedLoc(vacancy2)-FixedLoc(vacancy1);
	PathData.Path.PutInBox(disp1);
	PathData.Path.PutInBox(disp2);
	double dist=sqrt(dot(disp1,disp1));
	if (dist<Grid.End){
	  int index=Grid.ReverseMap(dist);
	  Histogram(index)++;
	}
	bool disp1Found=false;
	bool disp2Found=false;
    
	for (int counter3=0;counter3<FixedLoc.size();counter3++){
	  dVec dispCompare=FixedLoc(counter3)-FixedLoc(0);
	  PathData.Path.PutInBox(dispCompare);
	  if (dot(disp1-dispCompare,disp1-dispCompare)<=1e-5){
	    disp1Found=true;
	    HistogramDisp(counter3)=HistogramDisp(counter3)+1.0;
	  }
	  if (dot(disp2-dispCompare,disp2-dispCompare)<=1e-5){
	  HistogramDisp(counter3)=HistogramDisp(counter3)+1.0;
	  disp2Found=true;
	  }
	  
	}
// 	if (!disp1Found)
// 	  cerr<<"Disp 1 not found: "<<disp1<<" "<<disp2<<" "
// 	      <<vacancy1<<" "<<vacancy2
// 	      <<PathData.Path(slice,vacancy1)<<" "
// 	    <<PathData.Path(slice,vacancy2)<<endl;
// 	if (!disp2Found)
// 	  cerr<<"Disp 2 not found: "<<disp2<<" "<<disp1<<" "
// 	      <<vacancy1<<" "<<vacancy2
// 	      <<PathData.Path(slice,vacancy1)<<" "
// 	      <<PathData.Path(slice,vacancy2)<<endl;
	
	
      }
      else {
	int vacancy1=Perm(0)-1;	
	TempVacancyLoc(slice)=vacancy1;
      }
    }  
    if (numEmptySites==1){
      int halfSlice=PathData.Path.TotalNumSlices/2;
      for (int slice=0;slice<PathData.Path.NumTimeSlices()/2;slice++){
	int i=DispFromASite(TempVacancyLoc(slice),TempVacancyLoc(slice+halfSlice));
	if (i!=-1)
	  HistogramDisp(i)=HistogramDisp(i)+1.0;
	i=DispFromASite(TempVacancyLoc(slice+halfSlice),TempVacancyLoc(slice));
	if (i!=-1)
	  HistogramDisp(i)=HistogramDisp(i)+1.0;
      }
    }
    NumSamples++;

    
}

void VacancyLocClass::WriteBlock()
{
  int nSlices=PathData.Path.TotalNumSlices;
  double norm = 1.0/((double)NumSamples*(double)nSlices);
  Array<double,1> histogramWrite(Histogram.size());
  for (int counter=0;counter<histogramWrite.size();counter++)
    histogramWrite(counter)=(double)Histogram(counter)*norm;
  HistogramVar.Write(histogramWrite);
  HistogramDisp=HistogramDisp*norm;
  HistogramDispVar.Write(HistogramDisp);
  HistogramVar.Flush();
  NumSamples = 0;
  Histogram=0;
  HistogramDisp=0.0;
}


void VacancyLocClass::Read(IOSectionClass &in)
{  

  int numFixedPoints;
  ObservableClass::Read(in);
  assert(in.ReadVar("NumFixedPoints",numFixedPoints));
  FixedLoc.resize(numFixedPoints);
  Histogram.resize(numFixedPoints);
  HistogramDisp.resize(numFixedPoints);
  Histogram=0;
  HistogramDisp=0.0;
  Array<double,2> positions;
  assert(in.ReadVar("LocationsToCompare",positions));
  assert(in.ReadVar("Frequency",Freq));
  //  assert(in.ReadVar("dumpFreq",DumpFreq));

  ///Verify you used the right number of points to compare against
  assert(positions.extent(0)==FixedLoc.size());
  assert(positions.extent(1)==NDIM);

  dVec pos;
  for (int loc=0;loc<FixedLoc.size(); loc++){
    for (int dim=0; dim<NDIM; dim++)
      pos(dim) = positions(loc,dim);
    FixedLoc(loc) = pos;
  }      
  

  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }

  Array<int,1> toDivide(Histogram.size());
  toDivide=0;
  for (int counter=0;counter<FixedLoc.size();counter++){
    for (int counter2=0;counter2<FixedLoc.size();counter2++){
      dVec disp=FixedLoc(counter)-FixedLoc(counter2);
      PathData.Path.PutInBox(disp);
      double dist=sqrt(dot(disp,disp));
      if (dist<Grid.End){
	int index=Grid.ReverseMap(dist);
	toDivide(index)++;
      }
    }
  }
  IOSection.WriteVar("Multiplicity",toDivide);


  TempVacancyLoc.resize(FixedLoc.size());

  //setting up the table that maps the displacement from the A site
  //Makes the assumption that FixedLoc(0) is part of the A site
  DispFromASite.resize(FixedLoc.size(),FixedLoc.size());
  DispFromASite=-1;
  for (int loc1=0;loc1<FixedLoc.size();loc1++){
    for (int loc2=0;loc2<FixedLoc.size();loc2++){
      dVec disp1=FixedLoc(loc2)-FixedLoc(loc1);
      PathData.Path.PutInBox(disp1);
      for (int counter3=0;counter3<FixedLoc.size();counter3++){
	dVec dispCompare=FixedLoc(counter3)-FixedLoc(0);
	PathData.Path.PutInBox(dispCompare);
	if (dot(disp1-dispCompare,disp1-dispCompare)<=1e-5){
	  DispFromASite(loc1,loc2)=counter3;
	}
      }
    }
  }
}



// // Fix to include final link between link M and 0
// void VacancyLocClass::Accumulate()
// {
//   cerr<<"Into vacancy accumulate"<<endl;
//   dVec displaceAmount;
//   double distanceAmount;
//   SiteEmptyAtSomeTimeSlice=0;
//   for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++){
//     TempLoc=0;
//     for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
//       double closestAmount=5*dot(PathData.Path.GetBox(),PathData.Path.GetBox());
//       int closestLoc=0;
//       for (int counter=0;counter<FixedLoc.size();counter++){
// 	displaceAmount=PathData.Path(slice,ptcl)-FixedLoc(counter);
// 	PathData.Path.PutInBox(displaceAmount);
// 	distanceAmount=dot(displaceAmount,displaceAmount);
// 	if (distanceAmount<closestAmount){
// 	  closestAmount=distanceAmount;
// 	  closestLoc=counter;
// 	}
//       }
//       //Commented out Mar 1      Loc(closestLoc)++;
//       TempLoc(closestLoc)++;
//       R2Dist+=closestAmount;
//     }
    
//     for (int counter=0;counter<TempLoc.size();counter++)
//       if (TempLoc(counter)==0 &&  NeighborsNotDoublyOccupied(counter))
// 	SiteEmptyAtSomeTimeSlice(counter)=1;
//     for (int counter=0;counter<TempLoc.size();counter++)
//       if (TempLoc(counter)==0 && NeighborsNotDoublyOccupied(counter)){
// 	Loc(counter)++;
// 	//	cerr<<"My temp Loc is "<<TempLoc<<endl;
//       }
//     //    bool doPair=true;
//     //    for (int counter=0;counter<TempLoc.size();counter++)
//     //      if (TempLoc(counter)>1){
//     //	doPair=false;
//     //      }
//     //    doPair=false;  //HACK! HACK! HACK!
//     //    if (doPair){
//     //      for (int counter=0;counter<TempLoc.size();counter++)
//     //	for (int counter2=0;counter2<TempLoc.size();counter2++)
//     //	  if (TempLoc(counter)==0 && TempLoc(counter2)==0){
//     //	    dVec disp=FixedLoc(counter)-FixedLoc(counter2);
//     //	    PathData.Path.PutInBox(disp);
//     //	    double dist=sqrt(dot(disp,disp));
//     //	    if (dist<Grid.End){
//     //	      int index=Grid.ReverseMap(dist);
//     //	      Histogram(index)++;
//     //	    }
//     //	  }
//     //    }
//     //  }
  

//     cerr<<"Pre-playing with stuff"<<FixedLoc.size()<<endl;    
//     Array<double,2> DistTable(1,1,ColumnMajorArray<2>());
//     DistTable.resize(FixedLoc.size(), FixedLoc.size());
//     DistTable=0.0;
//     cerr<<"a"<<endl;
//     for (int latticeSite=0;latticeSite<FixedLoc.size();latticeSite++){
//       for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
// 	//	cerr<<"Slice: "<<slice<<" "<<latticeSite<<" "<<ptcl<<endl;
// 	dVec disp;
// 	double dist2;
// 	//	PathData.Path.DistDisp(slice,latticeSite,ptcl,dist,disp);
// 	disp=PathData.Path(slice,ptcl)-FixedLoc(latticeSite);
// 	PathData.Path.PutInBox(disp);
// 	dist2=dot(disp,disp);
// 	//	cerr<<ptcl+2<<" "<<latticeSite<<endl;
// 	DistTable(latticeSite,ptcl+2)=dist2;
//       }
//     }
//     int n =FixedLoc.size();
//     //    *n=FixedLoc.size();
//     Array<int,1> Perm;
//     Perm.resize(FixedLoc.size());
//     cerr<<"Pre-fortran"<<endl;
//     F77_LSAPR (&n,DistTable.data(),Perm.data());
//     cerr<<"Post fortran"<<endl;
//     cerr<<"Open Sites: "<<Perm(0)-1<<" "<<Perm(1)-1<<endl;
//     if (PathData.Path.NumSpecies()>1){ //do vacancy-helium3 correlation
//       for (int counter=0;counter<TempLoc.size();counter++)
// 	if (TempLoc(counter)==0 && NeighborsNotDoublyOccupied(counter)){
// 	  dVec disp=FixedLoc(counter)-PathData.Path(slice,PathData.Path.Species(1).FirstPtcl);
// 	  PathData.Path.PutInBox(disp);
// 	  double dist=sqrt(dot(disp,disp));
// 	  if (dist<Grid.End){
// 	    int index=Grid.ReverseMap(dist);
// 	    Histogram(index)++;
// 	  }
// 	}
    

//     }
//     else { //do vacancy-vacancy correlation
//   ///TRying again 
//       cerr<<"Lattice Sites Other way: ";
//     for (int counter=0;counter<TempLoc.size();counter++)
//       for (int counter2=0;counter2<TempLoc.size();counter2++)
// 	if (TempLoc(counter)==0 && TempLoc(counter2)==0 &&
// 	    NeighborsNotDoublyOccupied(counter) && 
// 	    NeighborsNotDoublyOccupied(counter2)){
// 	  cerr<<counter<<" "<<counter2<<endl;
// 	  dVec disp=FixedLoc(counter)-FixedLoc(counter2);
// 	  PathData.Path.PutInBox(disp);
// 	  double dist=sqrt(dot(disp,disp));
// 	  if (dist<Grid.End){
// 	    int index=Grid.ReverseMap(dist);
// 	    Histogram(index)++;
// 	  }
// 	}   
//     }

//   }
//   for (int counter=0;counter<TempLoc.size();counter++)
//     if (SiteEmptyAtSomeTimeSlice(counter)==1)
//       NumEmptyLatticeSites=NumEmptyLatticeSites+1;
  


//   NumSamples++;
// }

// void VacancyLocClass::WriteBlock()
// {
//   double norm = 1.0/((double)NumSamples);
//   double R2Avg=R2Dist*norm;
//   //  cerr<<"About to write the double"<<endl;
//   //  R2Var.Write(R2Avg);
//   //  cerr<<"Written teh double"<<endl;

//   Array<int,1> locSum(Loc.size());
//   Array<double,1> locWrite(Loc.size());
//   locSum=0;
//   PathData.Path.Communicator.Sum(Loc,locSum);
//   for (int counter=0;counter<locSum.size();counter++){
//     locWrite(counter)=locSum(counter)*norm;
//   }
//   VacancyLocVar.Write(locWrite);
//   NumEmptyLatticeSitesVar.Write(NumEmptyLatticeSites/(double)NumSamples);
//   Array<double,1> histogramWrite(Histogram.size());
//   for (int counter=0;counter<histogramWrite.size();counter++)
//     histogramWrite(counter)=Histogram(counter)*norm;
//   HistogramVar.Write(histogramWrite);
//   HistogramVar.Flush();
//   Loc=0;
//   TempLoc=0;
//   NumSamples = 0;
//   R2Dist=0.0;
//   Histogram=0;
//   NumEmptyLatticeSites=0.0;
//   //  cerr<<"I'm done with that"<<endl;
// }

// void VacancyLocClass::Read(IOSectionClass &in)
// {  


//   int numFixedPoints;
//   ObservableClass::Read(in);
//   assert(in.ReadVar("NumFixedPoints",numFixedPoints));
//   Loc.resize(numFixedPoints);
//   Loc=0;
//   TempLoc.resize(numFixedPoints);
//   SiteEmptyAtSomeTimeSlice.resize(numFixedPoints);
//   TempLoc=0;
//   FixedLoc.resize(numFixedPoints);
//   Neighbors.resize(numFixedPoints);
//   Array<double,2> positions;
//   assert(in.ReadVar("LocationsToCompare",positions));

//   ///Verify you used the right number of points to compare against
//   assert(positions.extent(0)==Loc.size());
//   assert(positions.extent(1)==NDIM);
//   dVec pos;

//   for (int loc=0;loc<FixedLoc.size(); loc++){
//     for (int dim=0; dim<NDIM; dim++)
//       pos(dim) = positions(loc,dim);
//     FixedLoc(loc) = pos;
//   }      
  

//   if (PathData.Path.Communicator.MyProc()==0){
//     WriteInfo();
//     IOSection.WriteVar("Type","Scalar");
//   }
//   Array<int,1> toDivide(Histogram.size());
//   for (int counter=0;counter<FixedLoc.size();counter++){
//     for (int counter2=0;counter2<FixedLoc.size();counter2++){
//       dVec disp=FixedLoc(counter)-FixedLoc(counter2);
//       PathData.Path.PutInBox(disp);
//       double dist=sqrt(dot(disp,disp));
//       if (dist<Grid.End){
// 	int index=Grid.ReverseMap(dist);
// 	toDivide(index)++;
//       }
//     }
//   }
//   IOSection.WriteVar("Multiplicity",toDivide);
//   TabulateNearbySites();
//   PrintNearbySites();
  
  
// }



void VacancyLocClass::WriteInfo()
{


}
