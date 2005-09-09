#ifndef IO_H
#define IO_H

#include "IOHDF5.h"
#include "IOASCII.h"

namespace IO {

  template<typename T> bool 
  IOTreeClass::WriteVar (string name, T val)
  {


  }
  template<typename T, int RANK> bool 
  IOTreeClass::WriteVar (string name, Array<T,RANK> &val)
  {


  }

}

#endif
