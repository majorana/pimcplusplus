#include "CoupledPermuteStageClass.h"


double CoupledPermuteStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}


void CoupledPermuteStageClass::InitBlock()
{


}


void CoupledPermuteStageClass::Read (IOSectionClass &in)
{
  Array<double,1> gamma;
  assert (in.ReadVar("Gamma", gamma));
  assert (gamma.size()==4);
  for (int i=0; i<4; i++) {
    Table1.Gamma[i] = gamma(i);
    Table2.Gamma[i] = gamma(i);
  }
}

bool CoupledPermuteStageClass::Attempt (int &slice1, int &slice2, 
				      Array<int,1> &activeParticles, double &prevActionChange)
{

  if (activeParticles(0)==-1){
  //   int sendProc;
//     int recvProc;
//     int numProcs=PathData.InterComm.NumProcs();
//     int myProc=PathData.InterComm.MyProc();
//     sendProc=(myProc+1) % numProcs;
//     recvProc=((myProc-1) + numProcs) % numProcs;
    //    Array<int,1> dummy(1);
    //    dummy(0)=0;
    //    Array<int,1> dummy2(1);
    //    dummy2(0)=0;
    //    PathData.InterComm.SendReceive(sendProc,dummy,recvProc,dummy2);

  
//   /////   cerr<<"Attempting the coupled permute stage "<<PathData.InterComm.MyProc()<<endl;
    Array<int,1> coupledWeight(1);
    Array<int,1> coupledWeightSend(1);
    coupledWeightSend(0)=PathData.Path.Weight;
    coupledWeight(0)=0;
    int sendProc;
    int recvProc;
    int numProcs=PathData.InterComm.NumProcs();
    int myProc=PathData.InterComm.MyProc();
    sendProc=(myProc+1) % numProcs;
    recvProc=((myProc-1) + numProcs) % numProcs;
//   Array<int,1> dummy(1);
//   dummy(0)=0;
//   Array<int,1> dummy2(1);
//   dummy2(0)=0;
    PathData.InterComm.SendReceive (sendProc, coupledWeightSend,
				    recvProc, coupledWeight);
//   PathData.InterComm.SendReceive(sendProc,dummy,recvProc,dummy2);
//   //  //  cerr<<"I am "<<myProc<<" Send is "<<sendProc<<" Receive is "<<recvProc<<endl;
//   //  sleep(5);
//   //  cerr<<"My weight is "<<coupledWeightSend(0)<<" and my friends weight is "<<coupledWeight(0)<<endl;
  
    Array<int,1> doEven;
    Array<int,1> dummy;
    doEven.resize(1);
    doEven=0;
    dummy.resize(1);
    if (coupledWeightSend(0)*coupledWeight(0)==1){
// //     //    cerr<<"Currently the same!"<<endl;
      if (0==myProc){
	double myRand=PathData.Path.Random.Local();
	if (myRand>0.5)
	  doEven=1;
	else 
	  doEven=0;
// //       //      PathData.InterComm.Send(sendProc,doEven);
	PathData.InterComm.SendReceive(sendProc,doEven,recvProc,dummy);
	if (doEven(0)==1){
	  Forw->OnlyOdd=false;
	  Rev->OnlyOdd=false;
	  Forw->OnlyEven=true;
	  Rev->OnlyEven=true;
	}
	else {
	  Forw->OnlyOdd=true;
	  Rev->OnlyOdd=true;
	  Forw->OnlyEven=false;
	  Rev->OnlyEven=false;
	}
      }
      if (1==myProc){
// ///      //      PathData.InterComm.Receive(sendProc,doEven);
	PathData.InterComm.SendReceive(sendProc,dummy,recvProc,doEven);
	if (doEven(0)==0){
	  Forw->OnlyOdd=false;
	  Rev->OnlyOdd=false;
	  Forw->OnlyEven=true;
	  Rev->OnlyEven=true;
	}
	else {
	  Forw->OnlyOdd=true;
	  Rev->OnlyOdd=true;
	  Forw->OnlyEven=false;
	  Rev->OnlyEven=false;
	}
      }
    }
    else {
// //     //    cerr<<"Currently different"<<endl;
      if (0==myProc){
	double myRand=PathData.Path.Random.Local();
	if (myRand>0.5)
	  doEven=1;
	else 
	  doEven=0;
// //       //      PathData.InterComm.Send(sendProc,doEven);
	PathData.InterComm.SendReceive(sendProc,doEven,recvProc,dummy);
	if (doEven(0)==0){
	  Forw->OnlyOdd=false;
	  Rev->OnlyOdd=false;
	  Forw->OnlyEven=true;
	  Rev->OnlyEven=true;
	}
	else {
	  Forw->OnlyOdd=true;
	  Rev->OnlyOdd=true;
	  Forw->OnlyEven=false;
	  Rev->OnlyEven=false;
	}
      }
      if (1==myProc){
// //       //      PathData.InterComm.Receive(sendProc,doEven);
	PathData.InterComm.SendReceive(sendProc,dummy,recvProc,doEven);
	if (doEven(0)==0){
	  Forw->OnlyOdd=false;
	  Rev->OnlyOdd=false;
	  Forw->OnlyEven=true;
	  Rev->OnlyEven=true;
	}
	else {
	  Forw->OnlyOdd=true;
	  Rev->OnlyOdd=true;
	  Forw->OnlyEven=false;
	  Rev->OnlyEven=false;
	}
      }
    }
  




    Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);
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
    return (acceptProb > psi);
  }
  else
    return true;    
}
