#ifndef ACTIONS_CLASS_H
#define ACTIONS_CLASS_H
#include "ShortRangeClass.h"
#include "LongRangeClass.h"
#include "LongRangeRPAClass.h"
#include "ShortRangePotClass.h"
#include "LongRangePotClass.h"
#include "KineticClass.h"
#include "NodalActionClass.h"

class ActionsClass
{
private:
  Array<PairActionFitClass*,1> PairArray;
  Array<PairActionFitClass*,2> PairMatrix;
  PathDataClass &PathData;
  int MaxLevels; //is this the right place for this?
  bool UseRPA;
public:
  // Actions
  KineticClass Kinetic;
  ShortRangeClass ShortRange;
  LongRangeClass LongRange;
  LongRangeRPAClass LongRangeRPA;
  Array<ActionBaseClass *,1> NodalActions;
  //NodalClass Nodal;
  //DiagonalClass Diagonal;
  //ImportanceSampleClass ImportanceSample;

  // Potentials
  ShortRangePotClass ShortRangePot;
  LongRangePotClass  LongRangePot;
  

  void Read(IOSectionClass &in);
  ActionsClass(PathDataClass &pathData) : 
    ShortRange(pathData,PairMatrix),
    ShortRangePot(pathData, PairMatrix),
    LongRange(pathData,PairMatrix,PairArray), 
    LongRangeRPA(pathData, PairMatrix, PairArray),
    LongRangePot(pathData, PairMatrix),
    Kinetic(pathData),
    PathData(pathData)
  {
    ///Do nothing for now
  }





};

#endif
