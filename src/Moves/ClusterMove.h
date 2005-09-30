#ifndef CLUSTER_MOVE_H
#define CLUSTER_MOVE_H

#include "MoveBase.h"

class LocalFlip : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  void MakeMove();
  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="LocalFlip");
      assert(moveInput.ReadVar("name",Name));

    }
  void AssignPtcl(int mol,Array<int,1>& activeParticles);
  double MolPairAction(int slice,int m,int n);
 // void RotateMol(int slice,int mol,dVec axis,double phi);
  void RotateMol(int slice,int mol,dVec Q);
  dVec ArbitraryRotate(dVec axis,dVec coord, double phi);
  void Strip(dVec R, dVec u,dVec &aligned, dVec &perp);
  void IntegrityCheck(int slice, Array<int,1> activeParticles);
  LocalFlip(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};

class GlobalFlip : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  void MakeMove();
  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="GlobalFlip");
      assert(moveInput.ReadVar("name",Name));

    }
  void AssignPtcl(int mol,Array<int,1>& activeParticles);
  double MolPairAction(int slice,int m,int n);
  void RotateMol(int slice,int mol,dVec Q);
  GlobalFlip(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(5);
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};

#endif
