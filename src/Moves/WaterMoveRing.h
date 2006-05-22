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

#ifndef WATER_MOVE_RING_H
#define WATER_MOVE_RING_H

#include "MoveBase.h"


class WaterRotateRing : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  void MakeMove();
  void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="WaterRotateRing");
      assert(moveInput.ReadVar("name",Name));

    }
  void  Molecule2Atoms(int moleculeNum);
  Array<int,1> coord_loc;
  dVec Rotate(dVec coord,int u1,int u2,int u3,double theta);
  WaterRotateRing(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
      /* do nothing for now */
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};


class WaterTranslateRing : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  void MakeMove();
  void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="WaterTranslateRing");
      assert(moveInput.ReadVar("name",Name));
      

    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };

  Array<int,1> coord_loc;
  void Molecule2Atoms(int moleculeNum);
  dVec Translate(double epsilon);
  WaterTranslateRing(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
      /* do nothing for now */
    }
  


};



#endif
