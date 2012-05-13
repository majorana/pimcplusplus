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

#ifndef GLOBAL_JOSEPHSON_MOVE_H
#define GLOBAL_JOSPEHSON_MOVE_H

#include "MoveBase.h"

class GlobalJosephsonMove : public ParticleMoveClass
{
 public:
  int nMax;
  void MakeMove();
  void MakeMoveSlow();
  double g(int j);
  double Ec;
  void BuildA();
  Array<double,1> A;
  void Read(IOSectionClass &moveInput);
  GlobalJosephsonMove(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection)
    {
      Ec=1;
      ActiveParticles.resize(1);
      nMax=5;
    }
  void WriteRatio()
    {
      //do nothing for now
    };
};


#endif
