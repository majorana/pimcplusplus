#ifndef PIMC_CLASS_H
#define PIMC_CLASS_H


#include "PathDataClass.h"
#include "Moves/MoveClass.h"
#include "Observables/ObservableClass.h"
#include "LoopClass.h"
#include "RunInfoClass.h"

class PIMCClass 
{

public:
  std::list<MoveClass*> Moves;
  std::list<ObservableClass*> Observables;
  void ReadMoves(IOSectionClass &in);
  void ReadObservables(IOSectionClass &in);
  void ReadAlgorithm(IOSectionClass &in);
  void WriteSystemInfo();
  string OutFileName;
  IOSectionClass OutFile;
  LoopClass Algorithm;
  RunInfoClass RunInfo;
public:
  PathDataClass PathData;
  void Read(IOSectionClass &in);
  void Run();
  PIMCClass() : Algorithm(*this, OutFile, Moves, Observables)
  {
    RunInfo.ProgramName="pimc++";
  }
};


#endif
