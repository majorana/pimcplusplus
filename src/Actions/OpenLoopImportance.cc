#include "OpenLoopImportance.h"
#include "../PathDataClass.h"




///This has to be called after pathdata knows how many
///particles it has
void OpenLoopImportanceClass::Read(IOSectionClass& in)
{
  string impChoiceString;
  if (PathData.Path.OpenPaths){
    assert(in.OpenSection("OpenLoop"));
    assert(in.ReadVar("ImportanceSample",impChoiceString));
    if (impChoiceString=="None")
      ImpChoice=NOIMP;
    else if (impChoiceString=="Distance")
      ImpChoice=DISTIMP;
    else if (impChoiceString=="Displacement")
      ImpChoice=DISPXIMP;
    else {
      cerr<<"You have given an invalid choice"<<endl;
      assert(1==2);
    }
    in.CloseSection();
  }
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
  //  int myProc=(PathData.InterComm.MyProc() % 16);
  int myProc=PathData.GetCloneNum();
  double shift=(myProc%16)+0.5;
  shift=10.5;
  if (ImpChoice==NOIMP)  
    return 0.0;
  else if (ImpChoice==DISTIMP)
    return -log(exp(-(dist-shift)*(dist-shift)));
  else if (ImpChoice==DISPXIMP)
    return -log(exp(-(disp(0)-shift)*(disp(0)-shift)));
  else {
    cerr<<"You haven't give a valid choice!"<<endl;
    assert(1==2);
  }

  //  cerr<<"MY shift is "<<shift<<endl;
  //  return -log(0.01+(dist*dist)*(0.94*exp(-dist*dist)+0.06));

  /////  return -log(exp(-0.5*(dist-3.5)*(dist-3.5))); //(dist*dist));

}



double OpenLoopImportanceClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  return 0.0;
}
