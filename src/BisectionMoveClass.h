#ifndef BISECTIONMOVE_CLASS_H
#define BISECTIONMOVE_CLASS_H

#include "Common.h"
#include "MoveClass.h";
#include "ActionClass.h"

class BisectionMoveClass : public MoveClass
{
private:

  int NumOfParticlesToMove;
  int ParticleTypeToMove;
  int StartTimeSlice;
  int NumLevels;
  
  void makeMove();
 public:
  BisectionMoveClass();
};







       
  














#endif
