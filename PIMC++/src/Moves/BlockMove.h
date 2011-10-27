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

#ifndef BLOCK_MOVE_H
#define BLOCK_MOVE_H

#include "MoveBase.h"
#include "PermuteTableClass.h"
#include "PermuteTableOnClass.h"
#include "BisectionClass.h"

class CycleBlockMoveClass  : public MoveClass
{
private:
  BisectionClass Bisection;
  PermuteTableClass Table1, Table2;
  PermuteTableClass *Forw, *Rev;
  int NumAccepted;
  int NumMoves;
public:
  int StepsPerBlock;
  int SpeciesNum;
  int NumLevels;
  void MakeMove();
  ///All moves ought to be able to read
  void Read(IOSectionClass &input);
  double AcceptanceRatio();

  CycleBlockMoveClass(PathDataClass &myPathData, IOSectionClass outSection ) : 
    MoveClass(myPathData, outSection), Bisection(myPathData),
    Table1(myPathData), Table2(myPathData)
  { 
    Forw = &Table1;
    Rev  = &Table2;
    NumAccepted=0;
    NumMoves=0;
  }
}
;



#endif
