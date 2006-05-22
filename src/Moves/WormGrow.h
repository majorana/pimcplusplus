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

#ifndef WORM_GROW_MOVE_H
#define WORM_GROW_MOVE_H

#include "WormStage.h"
#include "MultiStage.h"
#include "../PathDataClass.h"



/// This is the displace move, which attempts to pick up a whole
/// single-particle path and displace it by a random amount.
class WormGrowMoveClass : public MultiStageClass
{
private:
  /// This is the standard distribution of the displacement gaussian
  double Sigma;
  WormStageClass WormGrowStage;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  int MaxGrowth;
  WormGrowMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection),
    WormGrowStage(pathData, outSection)
  {
    DumpFreq = 20;
  }
};


#endif
