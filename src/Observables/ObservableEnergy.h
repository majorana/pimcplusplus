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

#ifndef OBSERVABLE_ENERGY_H
#define OBSERVABLE_ENERGY_H

#include "ObservableBase.h"

class EnergyClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum;//, TotalActionSum, ExpTotalActionSum, TIP5PSum;

  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar;//, TotalActionVar, ExpTotalActionVar, TIP5PVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  EnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator),
      KineticVar("Kinetic",IOSection,myPathData.Path.Communicator), 
      dUShortVar("dUShort",IOSection,myPathData.Path.Communicator), 
      dULongVar ("dULong", IOSection,myPathData.Path.Communicator), 
      NodeVar   ("Node",   IOSection,myPathData.Path.Communicator), 
      VShortVar ("VShort",IOSection,myPathData.Path.Communicator), 
      VLongVar  ("VLong",IOSection,myPathData.Path.Communicator)
    // TotalActionVar ("TotalAction",IOSection,myPathData.Path.Communicator),
    // ExpTotalActionVar ("ExpTotalAction",IOSection,myPathData.Path.Communicator)
    // TIP5PVar  ("TIP5P",IOSection,myPathData.Path.Communicator)
  {
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
//     TotalActionSum = 0.0;
//     ExpTotalActionSum = 0.0;
    NumSamples = 0;
    TimesCalled=0;
    //    TIP5PSum = 0;
  }
};

class EnergySignClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum;

  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  EnergySignClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator),
      KineticVar("Kinetic",IOSection,myPathData.Path.Communicator), 
      dUShortVar("dUShort",IOSection,myPathData.Path.Communicator), 
      dULongVar ("dULong", IOSection,myPathData.Path.Communicator), 
      NodeVar   ("Node",   IOSection,myPathData.Path.Communicator), 
      VShortVar ("VShort",IOSection,myPathData.Path.Communicator), 
      VLongVar  ("VLong",IOSection,myPathData.Path.Communicator)
  {
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
