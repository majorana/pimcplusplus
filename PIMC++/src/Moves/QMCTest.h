#ifndef QMCTEST_MOVE_H
#define QMCTEST_MOVE_H

#include "MoveBase.h"

class QMCTestMove : public ParticleMoveClass
{
  public:
  int numAccepted,numMoves;
  void MakeMove();
  void AssignPtcl(int mol,Array<int,1>& activeParticles);
  dVec Translate(double epsilon);
  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("Type",typeCheck));
      assert(typeCheck=="QMCTestMove");
      assert(moveInput.ReadVar("Name",Name));
    }

  QMCTestMove(PathDataClass &myPathData,IOSectionClass outSection) : 
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

#endif
