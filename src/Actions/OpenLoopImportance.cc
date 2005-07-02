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
    cerr<<"The ImportanceSampling string was "<<impChoiceString<<endl;
    if (impChoiceString=="None")
      ImpChoice=NOIMP;
    else if (impChoiceString=="Distance")
      ImpChoice=DISTIMP;
    else if (impChoiceString=="Displacement")
      ImpChoice=DISPXIMP;
    else if (impChoiceString=="Tuned")
      ImpChoice=TUNEDFUNCTION;
    else {
      cerr<<"You have given an invalid choice"<<endl;
      assert(1==2);
    }
    if (ImpChoice!=NOIMP){
      string shiftTypeString;
      assert(in.ReadVar("ShiftType",shiftTypeString));
      if (shiftTypeString=="FixedShift")
	assert(in.ReadVar("ShiftAmount",Shift));
      else if (shiftTypeString=="ProcShift"){
	int myProc=PathData.GetCloneNum();
	Shift=(myProc%16)+0.5;
      }
      else {
	cerr<<"I don't know what shift you want me to use"<<endl;
	assert(1==2);
      }
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
  //  int myProc=PathData.GetCloneNum();
  //  double shift=(myProc%16)+0.5;
  //  shift=6.5;
  if (ImpChoice==NOIMP){ 
    //    cerr<<"I have chosen no importance function"<<endl;
    //    assert(1==3);
    //    return 0.0;
    return -log(1.0/(dist*dist+0.05)); 
  }
  else if (ImpChoice==DISTIMP){
    //    cerr<<"I have chosen a distance importance function"<<endl;
    return -log(exp(-(dist-Shift)*(dist-Shift)));
  }
  else if (ImpChoice==DISPXIMP){
    //    cerr<<"I have chosen a displacement importance function"<<endl;
    return -log(exp(-(disp(0)-Shift)*(disp(0)-Shift)));
  }
  else if (ImpChoice==TUNEDFUNCTION){
    return log(10.0)*(-0.24*dist*dist+0.007*dist*dist*dist+0.0006*dist*dist*dist*dist)+-log(1.0/(dist*dist+0.05));
  }
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
