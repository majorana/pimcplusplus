#include "JosephsonPathDump.h"



void JosephsonPathDumpClass::Accumulate()
{

}


void JosephsonPathDumpClass::WriteBlock()
{
  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices();
  int numTimeSlices=PathData.Path.NumTimeSlices();
  PathClass &Path = PathData.Path;
  for (int slice=slice1;slice<slice2;slice++)
    Phase(slice)=Path(slice,0)[0];
  PhaseVar.Write(Phase);  
}

void JosephsonPathDumpClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}

