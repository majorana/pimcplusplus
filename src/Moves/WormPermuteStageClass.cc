#include "PermuteStageClass.h"
#include "WormPermuteStageClass.h"

double WormPermuteStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}

void WormPermuteStageClass::Accept()
{
  NumAccepted(0)++;
  NumAttempted(0)++;
}

void WormPermuteStageClass::Reject()
{
  NumAttempted(0)++;
}


void WormPermuteStageClass::WriteRatio()
{
//   Array<int,1> numAttemptTotal(4);
//   Array<int,1> numAcceptTotal(4);
//   Array<double,1> ratioTotal(4);
//   int totalAttempts=0;
//   cerr<<"A"<<endl;
//   PathData.Path.Communicator.Sum(NumAttempted,numAttemptTotal);
//   cerr<<"A2"<<endl;
//   PathData.Path.Communicator.Sum(NumAccepted,numAcceptTotal);
//   for (int len=0;len<4;len++){
//     totalAttempts=totalAttempts+numAcceptTotal(len);
//     if (numAttemptTotal(len)!=0)

//       ratioTotal(len)=(double)numAcceptTotal(len)/((double)numAttemptTotal(len));
//     else
//       ratioTotal(len)=0.0;
//   }
//   cerr<<"B"<<endl;
//       ///divides by 2because accept gets called twice in acepting stages
//   for (int i=0;i<numAttemptTotal.size();i++)
//     numAttemptTotal(i)=(int)(numAttemptTotal(i)/2);
//   cerr<<"C"<<endl;
 
//   if (totalAttempts!=0){
//     AcceptanceRatioVar.Write(ratioTotal);
//     AcceptanceTotalVar.Write(numAttemptTotal);
//     NumAttempted=0;
//     NumAccepted=0;
//   }
//   cerr<<"out of my permute stage ratio"<<endl;

}

void WormPermuteStageClass::InitBlock()
{


}


void WormPermuteStageClass::Read (IOSectionClass &in)
{
      
}

bool WormPermuteStageClass::Attempt (int &slice1, int &slice2, 
				      Array<int,1> &activeParticles, double &prevActionChange)
{
  if (activeParticles(0)==-1){

    while (!(PathData.Path.OpenLink-(1<<NumLevels)>=1 && PathData.Path.OpenLink<=PathData.Path.NumTimeSlices()-1)){
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
    }

    PathData.MoveJoin(0);
    slice2=PathData.Path.OpenLink;
    slice1=PathData.Path.OpenLink-(1<<NumLevels);
    int permuteWith=PathData.Path.OpenPtcl;
    while (permuteWith==PathData.Path.OpenPtcl)
      permuteWith=PathData.Path.Random.LocalInt(3);
    activeParticles.resize(2);
    activeParticles(0)=PathData.Path.OpenPtcl;
    activeParticles(1)=permuteWith;
    dVec tempPos=PathData.Path(slice1,permuteWith);
    PathData.Path(slice1,permuteWith)=PathData.Path(slice1,PathData.Path.OpenPtcl);
    PathData.Path(slice1,PathData.Path.OpenPtcl)=tempPos;
    
    int oldPerm=PathData.Path.Permutation(PathData.Path.OpenPtcl);
    PathData.Path.Permutation(PathData.Path.OpenPtcl)=PathData.Path.Permutation(permuteWith);
    PathData.Path.Permutation(permuteWith)=oldPerm;

    dVec oldPos=PathData.Path(slice1,PathData.Path.OpenPtcl);
    dVec newPos;///was /10 instead of /40 for the free particles
    newPos(0)= oldPos[0]+
      (PathData.Path.Random.Local()-0.5)*(PathData.Path.GetBox()(0)/5.0);
    newPos(1)= oldPos[1]+
      (PathData.Path.Random.Local()-0.5)*(PathData.Path.GetBox()(1)/5.0);
    newPos(2)= oldPos[2]+
      (PathData.Path.Random.Local()-0.5)*(PathData.Path.GetBox()(2)/5.0);
    //    newPos(0)=oldPos[0];
    //    newPos(1)=oldPos[1];
    //    newPos(2)=oldPos[2];

    PathData.Path.SetPos(PathData.Path.OpenLink,PathData.Path.NumParticles(),newPos);
    prevActionChange=0;
    return true;

  }
  else{


    return true;    
  }
}
