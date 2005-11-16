#include "EmptyStage.h"
void EmptyStageClass::Read(IOSectionClass &in)
{
  LocalStageClass::Read(in);
}

void EmptyStageClass::WriteRatio()
{
  //  cerr<<"About to write my ratio"<<endl;
  Array<double,1> acceptRatio(2);
  acceptRatio(0)=(double)AcceptRatio(0)/EndAttempts;
  acceptRatio(1)=(double)AcceptRatio(1)/EndAttempts;
  AcceptRatioVar.Write(acceptRatio);
  AcceptRatioVar.Flush();
  //  cerr<<"done writing my ratio"<<endl;
}


void EmptyStageClass::Accept()
{
  EndAttempts++;
}

void EmptyStageClass::Reject()
{
  EndAttempts++;
}



double EmptyStageClass::Sample(int &slice1,int &slice2, 
	      Array<int,1> &activeParticles)
{
  return 1.0;
}
