#ifndef PIMC_CLASS_H
#define PIMC_CLASS_H

#include "PathDataClass.h"
#include "MoveClass.h"
#include "ObservableClass.h"
#include "WrapClass.h"
#include "PermuteTableClass.h"


class PIMCClass 
{
private:
  Array<MoveClass*,1> Moves;
  Array<ObservableClass* ,1> Observables;
  void ReadMoves(IOSectionClass &in);
  void ReadObservables(IOSectionClass &in);
  void ReadAlgorithm(IOSectionClass &in);
  string OutFileName;
  IOSectionClass OutFile;
  LoopClass Algorithm;
public:
  PermuteTableClass PermuteTable;
  PathDataClass PathData;
  void Read(IOSectionClass &in);
  void Run();
  PIMCClass() : Algorithm(&Moves, &Observables), PermuteTable(PathData)
  { /* Do nothing for now */ }
};


#endif
