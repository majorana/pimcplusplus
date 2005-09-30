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
    else if (impChoiceString=="ReTuned")
      ImpChoice=RETUNEDFUNCTION;
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
  //    if (PathData.Path.Communicator.MyProc()==1)
  //    cerr<<"Into importance action"<<endl;

  int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
  if (procWithRefSlice != PathData.Path.Communicator.MyProc()) {
    //    if (PathData.Path.Communicator.MyProc()==1)
    //      cerr<<"leaving imp action with 0"<<endl;

    return 0.0;
  }

  int openLink=PathData.Path.OpenLink;
  int openPtcl=PathData.Path.OpenPtcl;
  //  if (PathData.Path.Communicator.MyProc()==1){
  //    cerr<<"not returning with 0"<<endl;
  //    cerr<<"Proc with ref slice is "<<procWithRefSlice<<endl;
  //    cerr<<"The open link is "<<openLink<<endl;
  //    cerr<<"The open ptcl is "<<openPtcl<<endl;
  //    cerr<<"The values are "<<PathData.Path(openLink,PathData.Path.NumParticles())<<" and "<<PathData.Path(openLink,openPtcl)<<endl;
  //  }
  dVec disp;
  double dist;
  PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			 dist,disp); //This is distance between head and tail!
  //  int myProc=(PathData.InterComm.MyProc() % 16);
  //  int myProc=PathData.GetCloneNum();
  //  double shift=(myProc%16)+0.5;
  //  //  shift=6.5;
  //    if (PathData.Path.Communicator.MyProc()==1)
  //      cerr<<"about to start hte ifs"<<endl;

  if (ImpChoice==NOIMP){ 
    //    if (PathData.Path.Communicator.MyProc()==1)
    //      cerr<<"Into none action"<<endl;

    //cerr<<"I have chosen no importance function"<<endl;
    //    assert(1==3);
    //    return 0.0;
    return -log(1.0/(dist*dist+0.05)); 
  }
  else if (ImpChoice==DISTIMP){
    //    if (PathData.Path.Communicator.MyProc()==1)
    //      cerr<<"Into dist action"<<endl;
    //    cerr<<"I have chosen a distance importance function"<<endl;
    return -log(exp(-(dist-Shift)*(dist-Shift)));
  }
  else if (ImpChoice==DISPXIMP){
    //    cerr<<"I have chosen a displacement importance function"<<endl;
    //    if (PathData.Path.Communicator.MyProc()==1)
    //      cerr<<"Into dispximp action"<<endl;
    int numLinks=PathData.Path.NumTimeSlices()-1;
    disp=0.0;
    for (int slice=0;slice<numLinks;slice++) {
      int realSlice=(openLink+slice) % numLinks;
      int realSlicep1=(openLink+slice+1) % numLinks;
      dVec linkDisp;
      linkDisp=PathData.Path.Velocity(realSlice,realSlicep1,openPtcl);
      disp =disp+linkDisp;
    }
    double  dist=sqrt(dot(disp,disp));
    return -6*dist;
    //    return -log(exp(-(disp(0)-Shift)*(disp(0)-Shift)));
  }
  else if (ImpChoice==TUNEDFUNCTION){
    //    if (PathData.Path.Communicator.MyProc()==1)
    //      cerr<<"Into tuned action"<<endl;

    return log(10.0)*(-0.24*dist*dist+0.007*dist*dist*dist+0.0006*dist*dist*dist*dist); ///HACK! only for OnlyX run+-log(1.0/(dist*dist+0.05));
  }
//   else if (ImpChoice==RETUNEDFUNCTION){
//     //    if (PathData.Path.Communicator.MyProc()==1)
//     //      cerr<<"Into retuned action"<<endl;

//     if (dist>9){
//       //      if (PathData.Path.Communicator.MyProc()==1)
//       //	cerr<<"leaving returned action"<<endl;
//       return 99999999;
//    }
  else if (ImpChoice==RETUNEDFUNCTION){
    //    if (PathData.Path.Communicator.MyProc()==1)
    //      cerr<<"Into retuned action"<<endl;

    //    if (dist>9){
      //      if (PathData.Path.Communicator.MyProc()==1)
      //	cerr<<"leaving returned action"<<endl;
    //      return 99999999;
      //    }
    //    else {
      //      if (PathData.Path.Communicator.MyProc()==1)
      //	cerr<<"now leaving retuned  action"<<endl;
      return log(10.0)*(-0.29*dist*dist+0.0012*dist*dist*dist*dist);
      //    }
  }
  else {
    cerr<<"You haven't give a valid choice!"<<endl;
    assert(1==2);
  }

  //  cerr<<"MY shift is "<<shift<<endl;
  //  return -log(0.01+(dist*dist)*(0.94*exp(-dist*dist)+0.06));

  /////  return -log(exp(-0.5*(dist-3.5)*(dist-3.5))); //(dist*dist));
  //    if (PathData.Path.Communicator.MyProc()==1)
  //      cerr<<"leaving ina  very weird way"<<endl;

}



double OpenLoopImportanceClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  return 0.0;
}
