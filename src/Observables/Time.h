/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef MCTIME_H
#define MCTIME_H

#include "ObservableBase.h"
#include "../Moves/MoveBase.h"
#include <time.h>

///Currently only owrks in serial mode
class MCTimeClass : public ObservableClass
{

private:
  //Must be initialized
  Array<double,1> MoveTime;
  ObservableVecDouble1 MoveTimeVar;
  Array<double,1> ObservableTime;
  ObservableVecDouble1 ObservableTimeVar;
  ObservableDouble TotalTimeVar;
  int StartTime;
  double TotalTime;
  list<MoveClass*>  &Moves;
  list<ObservableClass*> &Observables;

public:
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  MCTimeClass(PathDataClass &myPathData, IOSectionClass &ioSection,
	      list<MoveClass*> &moves,
	      list<ObservableClass*> &observables)
    : ObservableClass(myPathData, ioSection) , 
      MoveTimeVar("MoveTime",IOSection,myPathData.Path.Communicator),
      ObservableTimeVar("ObservableTime",
			IOSection,myPathData.Path.Communicator),
      TotalTimeVar("TotalTime",
		   IOSection,myPathData.Path.Communicator),

      Moves(moves),
      Observables(observables)
  {
    TimesCalled=0;
  }
};


#endif 
