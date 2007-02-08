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

#ifndef EVENT_CLASS_H
#define EVENT_CLASS_H

#include <Common/IO/IO.h>

using namespace IO;

class PathDataClass;
class PathClass;
#ifdef BUILD_DEV
  class PathClassDev;
#endif

/// EventClass is the parent of ObservableClass and MoveClass.  It is
/// used as a handle by the algorithm to call the different
/// observables and moves.
class EventClass
{
protected:
  /// Stores the number of times called
  int TimesCalled;
  
  /// This stores the location in the output file where this
  /// observable or move puts its output
  IOSectionClass IOSection;
  
  /// Stores a reference to PathData
  PathDataClass &PathData;
  /// And a reference to Path for convenience
#ifdef BUILD_DEV
  PathClassDev &Path;
#else
  PathClass &Path;
#endif

public:
  /// Stores the time spent doing the event in seconds.
  double TimeSpent;

  /// Stores the name that is used in the algorithm to refer to this
  /// event
  string Name;

  virtual void DoEvent()=0;
  virtual void Read(IOSectionClass& IO)=0;
  EventClass(PathDataClass &pathData, IOSectionClass& out);
 };

#endif
