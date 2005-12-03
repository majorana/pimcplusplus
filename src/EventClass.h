#ifndef EVENT_CLASS_H
#define EVENT_CLASS_H

#include <Common/IO/IO.h>

using namespace IO;

class PathDataClass;

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
public:
  /// Stores the time spent doing the event in seconds.
  double TimeSpent;

  /// Stores the name that is used in the algorithm to refer to this
  /// event
  string Name;

  virtual void DoEvent()=0;
  virtual void Read(IOSectionClass& IO)=0;
  EventClass(PathDataClass &pathData, IOSectionClass& out) : 
    PathData(pathData), IOSection(out), TimeSpent(0.0), TimesCalled(0)
  {
    // do nothing else for now
  }
 };

#endif
