#include "PermuteStageClass.h"


double TablePermuteStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}

void TablePermuteStageClass::Accept()
{
  int myLen=Forw->CurrentCycle.Length;
  //  assert(myLen<=4);
  //  assert(myLen>0);
  NumAccepted(myLen-1)++;
  NumAttempted(myLen-1)++;
  PermuteTableClass* temp;
  if (myLen!=1){
    temp=Forw;
    Forw=Rev;
    Rev=temp;
  }
  NeedToRebuildTable=false;
}

void TablePermuteStageClass::Reject()
{
  int myLen=Forw->CurrentCycle.Length;
  //  assert(myLen<=4);
  //  assert(myLen>0);
  NumAttempted(myLen-1)++;
  NeedToRebuildTable=false;

}


void TablePermuteStageClass::WriteRatio()
{
//   cerr<<"Writing my permtute stage ratio"<<PathData.Path.Communicator.MyProc()<<endl;
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

void TablePermuteStageClass::InitBlock()
{
  NeedToRebuildTable=true;

}


void TablePermuteStageClass::Read (IOSectionClass &in)
{
  Array<double,1> gamma;
  double epsilon;
  assert (in.ReadVar("Gamma", gamma));
  assert (in.ReadVar("epsilon",epsilon));
  assert (gamma.size()==4);
  for (int i=0; i<4; i++) {
    Table1.Gamma[i] = gamma(i);
    Table2.Gamma[i] = gamma(i);
  }
  Table1.epsilon=epsilon;
  Table2.epsilon=epsilon;
      
}

bool TablePermuteStageClass::Attempt (int &slice1, int &slice2, 
				      Array<int,1> &activeParticles, double &prevActionChange)
{
  //   // First, decide on the chunk of slices we're working on
// //   int slice1;
// //   int slice2;
// //   int numSlices=PathData.NumTimeSlices();

// //   if (!PathData.Path.OpenPaths){
// //     double xi=PathData.Path.Random.Local();
// //     slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
// //     slice1=slice2-(1<<NumLevels);
// //   }
// //   else {
// //     // First, decide on the chunk of slices we're working on
// //     do{
// //       double xi=PathData.Path.Random.Local();
// //       slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
// //       slice1=slice2-(1<<NumLevels);
// //       //      if (slice2>=PathData.Path.NumTimeSlices()-1){
// // 	//	cerr<<"ERROR!"<<slice2<<endl;
// //       //	slice2--;
// //       //      }
// //     } while ((slice1<=(int)(PathData.Path.OpenLink) && (int)(PathData.Path.OpenLink)<=slice2));
// //   }
// //   //  if (slice2>=PathData.Path.NumTimeSlices()-1){
// //   //    cerr<<"ERROR! ERROR! ERROR!";
// //   //  }
// //   //  cerr<<slice1<<" "<<slice2<<" "<<PathData.Path.NumTimeSlices()<<endl;
// //   PathData.MoveJoin(slice2);
  
// //   int step = 0;
// //   // Now, construct the Forward table

//  cerr<<"Starting TablePermuteStageClass"<<endl;
  //  sleep(10);
  if (activeParticles(0)==-1){
    if ((PathData.Path.OpenPaths && slice1<=PathData.Path.OpenLink && 
	PathData.Path.OpenLink<=slice2) || 
	(PathData.Path.OpenLink==PathData.Path.NumTimeSlices()-1 &&
	 (slice1==0 || slice2==0))){
      if (NeedToRebuildTable)
	Forw->ConstructCycleTable(SpeciesNum, slice1, slice2,
				  PathData.Path.OpenPtcl);
    }
    else {
      if (NeedToRebuildTable)
	Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);
    //    cerr<<"Time slices chosen"<<endl;
    //    sleep(10);
    }
    int NumPerms = 0;
    // Choose a permutation cycle
    double forwT = Forw->AttemptPermutation();
    activeParticles.resize (Forw->CurrentCycle.Length);
    activeParticles = Forw->CurrentParticles();
    double revT = Rev->CalcReverseProb(*Forw);
    double Tratio = forwT/revT;
    int len=Forw->CurrentCycle.Length;
    if (len % 2 ==0){
      PathData.Path.Weight=PathData.Path.Weight*-1;
    }
    double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
    double psi = PathData.Path.Random.Local();

    Array<int,1> currentParticles=Forw->CurrentParticles();  
    double pi_ratio = exp(-actionChange+prevActionChange);
    double acceptProb = min(1.0, pi_ratio/Tratio);
    prevActionChange = actionChange;
    //    cerr<<"Ending TablePermuteStageClass"<<endl;
    //    sleep(10);

    return (acceptProb > psi);
  }
  else{
    //    cerr<<"Ending at an odd spot"<<endl;
    //    sleep(10);

    return true;    
  }
}
