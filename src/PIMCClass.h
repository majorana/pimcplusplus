#ifndef PIMC_CLASS_H
#define PIMC_CLASS_H

#include "PathDataClass.h"
#include "MoveClass.h"
#include "Observables/ObservableClass.h"
#include "Observables/ObservableEnergy.h"
#include "WrapClass.h"
#include "PermuteTableClass.h"
#include "RunInfoClass.h"


class PIMCClass 
{

public:
  Array<MoveClass*,1> Moves;
  Array<ObservableClass* ,1> Observables;
  void ReadMoves(IOSectionClass &in);
  void ReadObservables(IOSectionClass &in);
  void ReadAlgorithm(IOSectionClass &in);
  void WriteSystemInfo();
  string OutFileName;
  IOSectionClass OutFile;
  LoopClass Algorithm;
  RunInfoClass RunInfo;
public:
  //  PermuteTableClass ForwPermuteTable, RevPermuteTable;
  PathDataClass PathData;
  void Read(IOSectionClass &in);
  void Run();
  PIMCClass() : Algorithm(&Moves, &Observables)
		//	ForwPermuteTable(PathData), RevPermuteTable(PathData)
    {

    RunInfo.ProgramName="pimc++";
    RunInfo.Version="0.1";
    }
};


#endif
