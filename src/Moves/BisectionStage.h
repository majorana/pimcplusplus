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

#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"

class BisectionStageClass : public LocalStageClass
{
public:
  void WriteRatio();
  //  void Read(IOSectionClass& IO);
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  Array<double,2> Correlated;
  bool UseCorrelatedSampling;
  void Accept();
  void Reject();
  int TotalLevels;
  void CalcShift(Array<int,1> &activeParticles,int slice);
  BisectionStageClass(PathDataClass &pathData, int level,
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection)
  { 
    //do nothing for now
    BisectionLevel = level;

  }
};

#endif
