#include "EndStageClass.h"

void EndStageClass::WriteRatio()
{
  Array<double,1> acceptRatio(2);
  acceptRatio(0)=(double)AcceptRatio(0)/EndAttempts;
  acceptRatio(1)=(double)AcceptRatio(1)/EndAttempts;
  AcceptRatioVar.Write(acceptRatio);
  AcceptRatioVar.Flush();
}


void EndStageClass::Accept()
{
  if (Open==HEAD)
    AcceptRatio(0)++;
  else 
    AcceptRatio(1)++;
  EndAttempts++;

}

void EndStageClass::Reject()
{
  EndAttempts++;

}


///Chooses the time slices and moves the join so that the join is in
///the correct place for that time slice.
void EndStageClass::ChooseTimeSlices(int &slice1,int &slice2)
{
  if (PathData.Path.Random.Local()>0.5)
    Open=HEAD;
  else
    Open=TAIL;
    

  if (Open==HEAD){
    slice1=(int)PathData.Path.OpenLink;
    slice2=(1<<NumLevels)+slice1;
    
    if (PathData.Path.Communicator.NumProcs()>1)
      return;
    ///Shift the time slices so that from the head there are
    ///2^numlevels slices available  
    while (slice2>=PathData.Path.NumTimeSlices() || slice1==0){
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      slice1=(int)PathData.Path.OpenLink;
      slice2=(1<<NumLevels)+slice1;
    }
  }
  else if (Open==TAIL){
    slice2=(int)PathData.Path.OpenLink;
    slice1=slice2-(1<<NumLevels);
    if (PathData.Path.Communicator.NumProcs()>1)
      return;

    ///Shift the time slices so that from the tail there are
    ///2^numlevels slices available  
    while (slice1<0 || slice2==PathData.Path.NumTimeSlices()-1){
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      slice2=(int)PathData.Path.OpenLink;
      slice1=slice2-(1<<NumLevels);
    }
  }
  else {
    cerr<<"ERROR! We don't know whether to do head or tail?"<<endl;
    assert(1==2);
  }
  //  PathData.MoveJoin(slice2);

}

///HACK! HACK! HACK!
double mysign(double num)
{
  if (num<0)
    return -1.0;
  else return 1.0;
}

double EndStageClass::Sample(int &slice1,int &slice2, 
	      Array<int,1> &activeParticles)
{
//   cerr<<"About to sleep"<<endl;
//   sleep(2000);
//   cerr<<"Done to sleep"<<endl;
  int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
  //  cerr<<"Here I am "<<endl;
  //If you don't own the open link then just choose some slices to work with
  if (procWithRefSlice != PathData.Path.Communicator.MyProc()) {
    //    cerr<<"In here"<<endl;
    int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    slice1 = PathData.Path.Random.LocalInt (numLeft);
    slice2 = slice1+sliceSep;
    return 1.0;
  }
//   sleep(200);
//   cerr<<"Before choose time slices"<<endl;

  ChooseTimeSlices(slice1,slice2);
  //  cerr<<"Ending choose time slice"<<endl;
  ///if you are running with time slice parallelization it doesn't
  ///shift the data so everyone doesn't have to know about it. It
  ///simply returns 0 forcing the move to be rejected prematurely early
  if ((slice1<=0 || slice2>=PathData.Path.NumTimeSlices()-1) &&
      PathData.Path.Communicator.NumProcs()>1)
    return 0.0;
  PathData.MoveJoin(slice2);
  ///The active particle should be a slice
  activeParticles.resize (1);
  activeParticles(0) = (int)PathData.Path.OpenPtcl;


  int changePtcl;
  if (Open==HEAD){
    //    cerr<<"Setting to head"<<endl;
    changePtcl=(int)PathData.Path.OpenPtcl;
  }
  else if (Open==TAIL){
    //    cerr<<"Setting to Tail"<<endl;
    changePtcl=PathData.Path.NumParticles();
  }
  else {
    cerr<<"I'm neither seeing a head or tail here! Warning! Warning!"<<endl;
    abort();
  }

  dVec oldPos=PathData.Path(PathData.Path.OpenLink,changePtcl);
  dVec newPos;///was /10 instead of /40 for the free particles
  newPos(0)= oldPos[0]+
    PathData.Path.Random.Local()*(PathData.Path.GetBox()(0)/40.0)*
    mysign(PathData.Path.Random.Local()-0.5);
  newPos(1)= oldPos[1]+
    PathData.Path.Random.Local()*(PathData.Path.GetBox()(1)/40.0)*
    mysign(PathData.Path.Random.Local()-0.5);
  newPos(2)=  oldPos[2]+
    PathData.Path.Random.Local()*(PathData.Path.GetBox()(2)/40.0)*
    mysign(PathData.Path.Random.Local()-0.5);
  if (fabs(newPos(0))>500 || fabs(newPos(1))>500 || fabs(newPos(2))>500){
    cerr<<"ERROR! ERROR! ERROR!"<<newPos(0)<<" "<<newPos(1)<<" "<<newPos(2)<<endl;
  } 
  PathData.Path.SetPos(PathData.Path.OpenLink,changePtcl,newPos);
  return 1.0;
}
