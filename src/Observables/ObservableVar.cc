#include "ObservableVar.h"


void ObservableVar::
Flush()
{
  if (Comm.MyProc()==0)
    Out.FlushFile();
}

