#ifndef BISECTIONMOVE_CLASS_H
#define BISECTIONMOVE_CLASS_H

#include "Common.h"
#include "MoveClass.h"
#include "ActionClass.h"

class BisectionMoveClass : public ParticleMoveClass
{

  
 public:
  int StartTimeSlice;
  int NumLevels;

  void makeMove();
  BisectionMoveClass();
};


class ShiftMove : public MoveClass
{
 public:
  int numTimeSlicesToShift;
  void makeMove();
};





       
  














#endif
