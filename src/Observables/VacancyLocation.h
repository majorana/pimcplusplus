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

#ifndef VACANCY_LOC_H
#define VACANCY_LOC_H

#include "ObservableBase.h"
#include <list>
#include <vector>

class VacancyLocClass : public ObservableClass
{

private:
  //Must be initialized
  Array<int,1> Histogram;
  Array<int,1> Loc;
  Array<int,1> TempLoc;
  Array<int,1> SiteEmptyAtSomeTimeSlice;
  LinearGrid Grid;
  ///This is the set of locations you should compare against to decide
  ///the location of the head and the tail
  Array<dVec,1> FixedLoc;
  vector<list<int> > Neighbors;
  ObservableVecDouble1 VacancyLocVar;
  ObservableVecDouble1 HistogramVar;
  ObservableDouble NumEmptyLatticeSitesVar;
  double NumEmptyLatticeSites;
  //  ObservableDouble R2Var;
  double R2Dist;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void TabulateNearbySites();
  void PrintNearbySites();
  bool NeighborsVacancyFree(int site);
  double VacancyDistance();
  bool NeighborsNotDoublyOccupied(int site);
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  VacancyLocClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      VacancyLocVar("y",IOSection,myPathData.Path.Communicator),//,
    HistogramVar("g(r)",IOSection,myPathData.Path.Communicator),
    NumEmptyLatticeSitesVar("EmptySites",IOSection,myPathData.Path.Communicator)
    //    R2Var("R2",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0; 
    TimesCalled=0;
    R2Dist=0.0;
    int numPoints=100;
    Grid.Init(0,PathData.Path.GetBox()[0],numPoints);
    Histogram.resize(numPoints);
    Histogram=0;
    NumEmptyLatticeSites=0.0;
  }
};


#endif 
