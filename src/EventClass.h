#ifndef EVENT_CLASS_H
#define EVENT_CLASS_H

#include "Common/IO/InputOutput.h"


class EventClass
{
 public:
  virtual void DoEvent()=0;
  virtual void Read(IOSectionClass& IO)=0;
  

};

#endif
