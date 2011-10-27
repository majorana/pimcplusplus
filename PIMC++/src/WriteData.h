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

#ifndef WRITE_DATA_CLASS_H
#define WRITE_DATA_CLASS_H

#include "EventClass.h"
class ObservableClass;
class MoveClass;


/// Note:  LoopClass is not allow to write any output to its out
/// IOsection. 
class WriteDataClass : public EventClass 
{
protected:
  list<MoveClass*> &Moves;
  list<ObservableClass*> &Observables;
 public:
  void DoEvent();
  void Read(IOSectionClass &IO);
  WriteDataClass(PathDataClass &pathData, IOSectionClass &out,
	    list<MoveClass*> &moves,
	    list<ObservableClass*> &observables) : 
    EventClass(pathData, out),
    Moves(moves), Observables(observables)
  {
    Name = "WriteData";
  }
};



#endif
