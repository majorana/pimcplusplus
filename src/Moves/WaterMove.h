#ifndef WATER_MOVE_H
#define WATER_MOVE_H

#include "MoveBase.h"
#include <time.h>

const int numSlots = 100;//(int)(resolution*sqrt(3.0)*2*M_PI) + 1;

class WaterRotate : public ParticleMoveClass
{
 public:
  double Theta;
  int numAccepted,numMoves;
  clock_t newtime,oldtime,elapsed,watch,watchstart,watchend;
  void MakeMove();
  void Read(IOSectionClass &moveInput);
  void  Molecule2Atoms(int moleculeNum);
  Array<int,1> coord_loc;
  dVec Rotate(dVec coord,int u1,int u2,int u3,double theta);
  WaterRotate(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
      /* do nothing for now */
      oldtime = clock();
      watch = 0;
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};

class WaterRareRotate : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  Array<int,4> AcceptLog;
  void MakeMove();
  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="WaterRareRotate");
      assert(moveInput.ReadVar("name",Name));

    }
  void  Molecule2Atoms(int moleculeNum);
  Array<int,1> coord_loc;
  dVec Rotate(dVec coord,int u1,int u2,int u3,double theta);
  dVec ArbitraryRotate(dVec axis,dVec coord, double phi);
  void Strip(dVec R, dVec u,dVec &aligned, dVec &perp);
  void IntegrityCheck(int slice, Array<int,1> activeParticles);
  int FindSlot(double angularmag,int flag);
  WaterRareRotate(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
      AcceptLog.resize(numSlots);
      AcceptLog = 0;
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};

class WaterTranslate : public ParticleMoveClass
{
 public:
  double Step;
  int numAccepted,numMoves;
  void MakeMove();
  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="WaterTranslate");
      assert(moveInput.ReadVar("name",Name));
      assert(moveInput.ReadVar("step",Step)); 

    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  inline void WriteRatio()
    {
      //do nothing for now
    };

  Array<int,1> coord_loc;
  void Molecule2Atoms(int moleculeNum);
  dVec Translate(double epsilon);
  WaterTranslate(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
      /* do nothing for now */
    }
  


};


#endif
