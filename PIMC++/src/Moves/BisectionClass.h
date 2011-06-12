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

#ifndef BISECTION_CLASS_H
#define BISECTION_CLASS_H

#include "../PathDataClass.h"


class BisectionClass 
{
public:
  Array<unsigned long,1> NumAccepted, NumRejected;

public:
  PathDataClass  &PathData;

  void AcceptanceRatios(Array<double,1> &ratios);
  bool Bisect(int startTimeSlice,int numLevels, Array<int,1> activeParticles);
  bool Bisect(int startTimeSlice,int numLevels, Array<int,1> activeParticles,
	      double permActionChange);

  ///This picks a new location in space for the particles in the
  ///particles Array at all of the time slices between startSlice and
  ///endSlice (at the appropriate skip for the level)
  double SamplePaths(int startSlice, int endSlice, Array<int,1> particles, int level);

  /// This calculates the sample probability for going from the state
  /// that is currently in the newMode of MirroredArrayClass to the
  /// state that is currently in oldMode of MirroredArrayClass 
  double LogSampleProb(int startSlice, int endSlice, Array<int,1> particles, int level);
  BisectionClass(PathDataClass &pathData) : PathData(pathData)
  { 
    NumAccepted.resize(pathData.Action.MaxLevels);
    NumRejected.resize(pathData.Action.MaxLevels);
    NumAccepted = 0;
    NumRejected = 0;
  }
/* Do nothing for now */
};

#endif
