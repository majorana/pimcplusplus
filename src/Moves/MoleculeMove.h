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

#ifndef MOLECULE_MOVE_H
#define MOLECULE_MOVE_H

#include "MoleculeMoveBase.h"
#include "../Moves/MoveUtils.h"

class MoleculeRotate : public MolMoveClass
{
  bool doAllSlices;
 public:
  double Theta;
  void Set(double setTheta);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  MoleculeRotate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			cerr << "MoleculeRotate constructor" << endl;
    }

  MoleculeRotate(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start)
    {
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);

  //void WriteRatio()
  //  {
  //  };
};


class BondStretch : public MolMoveClass
{
 public:
  double s;
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  BondStretch(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			cerr << "BondStretch constructor" << endl;
    }

  BondStretch(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start)
    {
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);

  //void WriteRatio()
  //  {
  //  };
};

class MoleculeTranslate : public MolMoveClass
{
  bool doAllSlices;
	int counter;
 	public:
  double Step;
  void Set(double setStep);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  //inline void WriteRatio()
  //  {
  //  };

  MoleculeTranslate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			counter = 0;
    }

  MoleculeTranslate(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start)
    {
    }
};

// individual intramolecular moves
class ParticleTranslate : public MolMoveClass
{
	int counter;
 	public:
  double Sigma;
  void Set(double setSigma);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  //inline void WriteRatio()
  //  {
  //  };

  ParticleTranslate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			counter = 0;
    }

  ParticleTranslate(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start)
    {
    }
};

class DimerMove : public MolMoveClass
{
 	public:
  double Step;
  void Set(double setStep);
  //void MakeMove();
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  DimerMove(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
    }

  DimerMove(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start)
    {
    }
};

class DummyEvaluate : public MolMoveClass
{
	int counter;
 	public:
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);
  //inline void WriteRatio()
  //  {
  //  };

  DummyEvaluate(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection)
    {
			counter = 0;
    }

  DummyEvaluate(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start)
    {
    }
};

class MoleculeMulti : public MolMoveClass
{
  MoleculeRotate Rotate;
  MoleculeTranslate Trans;
  BondStretch Stretch;

 public:
  double Theta;
  double Step;
  double s;
  double Sample(int &slice1,int &slice2, Array<int,1> &activeParticles);
  void Read(IOSectionClass &moveInput);

  MoleculeMulti(PathDataClass &myPathData,IOSectionClass outSection) : 
    MolMoveClass (myPathData,outSection),
    Rotate(myPathData, outSection), Trans(myPathData, outSection), Stretch(myPathData, outSection)
    {
			cerr << "MoleculeMulti constructor" << endl;
    }

  MoleculeMulti(PathDataClass &myPathData,IOSectionClass outSection, int numToRead, int start) : 
    MolMoveClass (myPathData,outSection, numToRead, start),
    Rotate(myPathData, outSection, numToRead, start), Trans(myPathData, outSection, numToRead, start), Stretch(myPathData, outSection, numToRead, start)
    {
    }
};

#endif
