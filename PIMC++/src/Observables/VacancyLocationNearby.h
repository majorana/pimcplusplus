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

#ifndef VACANCY_LOC_NEARBY_H
#define VACANCY_LOC_NEARBY_H

#include "ObservableBase.h"
#include <list>
#include <vector>

class VacancyLocNearbyClass : public ObservableClass
{

private:
  //Must be initialized
  Array<int,1> Histogram;
  Array<double,2> HistogramDisp;
  LinearGrid Grid;
  ///This is the set of locations you should compare against to decide
  ///the location of the head and the tail
  Array<dVec,1> FixedLoc;
  Array<int,2> DispFromASite;
  Array<int,1> TempVacancyLocNearby;
  ObservableVecDouble1 HistogramVar;
  ObservableVecDouble2 HistogramDispVar;
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
  void WriteBlockv2();
  void Read(IOSectionClass& in);
  VacancyLocNearbyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      HistogramDispVar("HistogramDisp",IOSection,myPathData.Path.Communicator),
      HistogramVar("g(r)",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0; 
    TimesCalled=0;
    int numPoints=100;
    Grid.Init(0,PathData.Path.GetBox()[0],numPoints);
  }
};


#endif 
