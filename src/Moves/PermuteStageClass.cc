#include "PermuteStageClass.h"


void PermuteStageClass::InitBlock()
{
  //do nothing for now
}
void PermuteStageClass::Read (IOSectionClass &in)
{
  //do nothing for now
}
void PermuteStageClass::Accept()
{
  //do nothing for now
}

void PermuteStageClass::Reject()
{
  //do nothing for now
}



// double TablePermuteStageClass::Sample (int &slice1, int &slice2,
// 		 Array<int,1> &activeParticles)
// {
  

// }





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
  Forw->ConstructCycleTable(SpeciesNum, slice1, slice2);
  
  int NumPerms = 0;
  
  // Choose a permutation cycle
  double forwT = Forw->AttemptPermutation();
  double revT = Rev->CalcReverseProb(*Forw);
  double Tratio = forwT/revT;
  int len=Forw->CurrentCycle.Length;
  double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
  double psi = PathData.Path.Random.Local();
  Array<int,1> currentParticles=Forw->CurrentParticles();  
  double pi_ratio = exp(-actionChange);
  double accept_prob = min(1.0, pi_ratio/Tratio);
  //    if (log(psi) < -(actionChange + log(Tratio))) {
  if  (accept_prob > psi) {
    bool acceptBisect = 
      Bisection.Bisect(slice1, NumLevels, currentParticles, actionChange);
    if (acceptBisect){
      PathData.AcceptMove(slice1,slice2,currentParticles);
      
      // We don't construct a new table for single-ptcl moves!
      if (Forw->CurrentCycle.Length!=1){
	NumPerms++;
	PermuteTableClass* tempPtr=Forw;
	Forw=Rev;
	Rev=tempPtr;
      }
      NumAccepted++;
    }
    else{
      PathData.RejectMove(slice1,slice2,currentParticles);
    }
  }
  else
    PathData.RejectMove(slice1, slice2,currentParticles);
  

 


}
