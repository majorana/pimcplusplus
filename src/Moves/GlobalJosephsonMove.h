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
