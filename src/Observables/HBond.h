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

#ifndef HBOND_H
#define HBOND_H

#include "ObservableBase.h"

class HbondClass : public ObservableClass
{

private:
  double TotalSum;

  ObservableDouble TotalVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  bool CheckPair(int slice, int obond, int ohome, int p);
  bool IsHBond(int slice, int OA, int OB);
  void Accumulate();
  void WriteBlock();
  void Initialize();
  void Read(IOSectionClass& in);
  HbondClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator)
  {
    Initialize();
  }
};
#endif 
