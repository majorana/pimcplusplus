#ifndef OBSERVABLE_BASE_H
#define OBSERVABLE_BASE_H


#include "../Common.h"
#include "../PathDataClass.h"
#include "../Common/IO/InputOutput.h"


/// This is the parent class for all observables.  It contains
/// a pointer to PathData.
class ObservableClass 
{
protected:
  /// The first time you write to an observable you have to do the
  /// write a little differently and you might need to write additional
  /// info like the description, etc.
  bool FirstTime;
  ///You can add more IOVar pointers to inhereted classes if necessary
  VarClass *IOVar;
public:
  /// A reference to the PathData I'm observing
  PathDataClass &PathData;
  /// Note: This is not a reference.  If it were, it could change
  /// behind our backs
  IOSectionClass IOSection;  
  string Name;
  string Description;
  /// Observe the state of the present path and add it to the
  /// running sum for averages.
  virtual void Accumulate() = 0;
  virtual void WriteBlock()=0;
  virtual void Read(IOSectionClass& IO);
  ////We don't actually ever call this so we really shouldn't allow
/// it's existence yet. 
  ///  virtual void ShiftData(int numTimeSlices) {;}
  virtual void WriteInfo();

  /// The constructor.  Sets PathData references and calls initialize.
  ObservableClass(PathDataClass &myPathData,IOSectionClass ioSection) 
    : PathData(myPathData), IOSection(ioSection)
  {
    FirstTime = true;
    Name="";
    Description="";
  }
};




#endif
