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

#ifndef AUTO_CORRELATION_H
#define AUTO_CORRELATION_H

#include "ObservableBase.h"

/// An autocorrelation observable.
class AutoCorrClass : public ObservableClass
{
  /// Stores number of counts in each bin
  Array<double,1> OneOneHistogram;
  Array<double,1> OneNetHistogram;
  Array<double,1> NetNetHistogram;
  ObservableVecDouble1 OneOneHistVar;
  ObservableVecDouble1 OneNetHistVar;
  ObservableVecDouble1 NetNetHistVar;
  ObservableDouble TotalNetMuVar;
	double TotalNetMu;
  /// Stores dipole moments
  Array<dVec,2> DipoleBin;
  Array<dVec,2> NetDipoleBin;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  //int Freq;
  int dumpFrequency;
  int now;
	string dipoleSpecies;
public:
	int NumSlots;
	int WaitToFill;
	int LastTotal;
  /// This grid defines the bins.  Bin 0 is bounded by 0 on the 
  /// bottom and grid(0) on the top.
  LinearGrid grid;
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void Print();
  void WriteBlock();
  void LocalWriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  int MapIndex(int slice, int molecule);
  void Advance(int& index,int limit);
  dVec MeasureDipole(int slice,int molecule);
  void CalcAutoCorr(int index, int t, int limit, double& SingleSingle, double& SingleNet, double& NetNet);
  int Locate(int i, int t, int limit);
  double CalcDotProd(dVec v1, dVec v2);
  dVec Rotate(dVec coord, double theta);
  AutoCorrClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    OneOneHistVar("OneOney", IOSection, myPathData.Path.Communicator),
    OneNetHistVar("OneNety", IOSection, myPathData.Path.Communicator),
    NetNetHistVar("NetNety", IOSection, myPathData.Path.Communicator),
    TotalNetMuVar  ("NetDipoleMoment",  IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
    Initialize();
  }

/*  AutoCorrClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData, ioSection) 
  { Initialize(); }
*/

};
#endif
