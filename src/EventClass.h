#ifndef EVENT_CLASS_H
#define EVENT_CLASS_H

#include <Common/IO/IO.h>

using namespace IO;

class EventClass
{
 public:
  virtual void DoEvent()=0;
  virtual void Read(IOSectionClass& IO)=0;
  

};

#endif
