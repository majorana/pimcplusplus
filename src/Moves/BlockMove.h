#ifndef BLOCK_MOVE_H
#define BLOCK_MOVE_H

#include "MoveBase.h"
#include "PermuteTableClass.h"
#include "PermuteTableOnClass.h"
#include "BisectionClass.h"

class CycleBlockMoveClass  : public MoveClass
{
private:
  BisectionClass Bisection;
  PermuteTableClass Table1, Table2;
  PermuteTableClass *Forw, *Rev;
  int NumAccepted;
  int NumMoves;
public:
  int StepsPerBlock;
  int SpeciesNum;
  int NumLevels;
  void MakeMove();
  ///All moves ought to be able to read
  void Read(IOSectionClass &input);
  double AcceptanceRatio();

  CycleBlockMoveClass(PathDataClass &myPathData, IOSectionClass outSection ) : 
    MoveClass(myPathData, outSection), Bisection(myPathData),
    Table1(myPathData), Table2(myPathData)
  { 
    Forw = &Table1;
    Rev  = &Table2;
    NumAccepted=0;
    NumMoves=0;
  }
}
;



#endif
