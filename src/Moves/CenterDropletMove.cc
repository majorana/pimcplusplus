#include "CenterDropletMove.h"


void CenterDropletClass::MakeMove()
{
  ofstream outfile;
  outfile.open("CenterOfMass.txt");
  PathData.MoveJoin(PathData.Path.NumTimeSlices()-1);
  dVec totalDispVec(0.0);
  
  int firstPtcl=PathData.Path.Species(Species).FirstPtcl;
  int lastPtcl=PathData.Path.Species(Species).LastPtcl;

  dVec chainCtr;
  //  cerr<<"My first particle is "<<firstPtcl<<endl;
  //  cerr<<"My last particle is "<<lastPtcl<<endl;
  for (int ptcl=firstPtcl;ptcl<=lastPtcl;ptcl++){
    dVec dispToPtcl=PathData.Path(0,ptcl);
    PathData.Path.PutInBox(dispToPtcl);
    chainCtr(0)=0.0;
    chainCtr(1)=0.0;
    chainCtr(2)=0.0;
    for (int slice=1;slice<PathData.Path.NumTimeSlices();slice++){
      dVec disp;
      disp=PathData.Path.Velocity(0,slice,ptcl);
      chainCtr+=disp;
    }
    //    PathData.Path.PutInBox(chainCtr);
    chainCtr = chainCtr/PathData.Path.NumTimeSlices();
    chainCtr=chainCtr+dispToPtcl;
    PathData.Path.PutInBox(chainCtr);
    //    cerr<<chainCtr[0]<<" "<<chainCtr[1]<<" "<<chainCtr[2]<<endl;
    outfile<<chainCtr[0]<<" "<<chainCtr[1]<<" "<<chainCtr[2]<<endl;
    totalDispVec+=chainCtr;
  }
  PathData.Path.PutInBox(totalDispVec);
  totalDispVec /= (lastPtcl-firstPtcl+1);
  //  cerr<<"My total disp vector is "<<totalDispVec<<endl;
  outfile.close();

 
  for (int ptcl=firstPtcl;ptcl<=lastPtcl;ptcl++)
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++)
      PathData.Path(slice,ptcl)-=totalDispVec;

  PathData.AcceptMove(0,PathData.Path.NumTimeSlices()-1,ActiveParticles);

}
