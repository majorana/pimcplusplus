#include "OpenLoopImportance.h"
#include "../PathDataClass.h"

///This has to be called after pathdata knows how many
///particles it has
void OpenLoopImportanceClass::Read(IOSectionClass& in)
{
}

OpenLoopImportanceClass::OpenLoopImportanceClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double OpenLoopImportanceClass::Action (int slice1, int slice2,
			     const Array<int,1> &changedParticles, int level)
{
  int openLink=PathData.Path.OpenLink;
  int openPtcl=PathData.Path.OpenPtcl;
  dVec disp;
  double dist;
  PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			 dist,disp); //This is distance between head and tail!
  int myProc=(PathData.InterComm.MyProc() % 16);
  double shift=myProc+0.5;
  //  return 0.0;
  return -log(exp(-(dist-shift)*(dist-shift)));

  //  cerr<<"MY shift is "<<shift<<endl;
  //  return -log(0.01+(dist*dist)*(0.94*exp(-dist*dist)+0.06));

  /////  return -log(exp(-0.5*(dist-3.5)*(dist-3.5))); //(dist*dist));

}



double OpenLoopImportanceClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  return 0.0;
}
