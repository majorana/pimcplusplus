#ifndef BISECTIONMOVE_CLASS_H
#define BISECTIONMOVE_CLASS_H

#include "Common.h"
#include "MoveClass.h"
#include "ActionClass.h"

class BisectionMoveClass : public ParticleMoveClass
{
private:



  int StartTimeSlice;
  int NumLevels;
  
  void makeMove();
 public:
  BisectionMoveClass();
};


class ShiftMove : public MoveClass
{
  int numTimeSlicesToShift;
  void makeMove();
  


}





       
  














#endif
