/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include <Common/MPI/Communication.h>
#include "MultiStage.h"
#include "sys/time.h"

void MultiStageClass::Read(IOSectionClass& in)
{
  cm2=0.0;
  ///do nothing for now
}

void MultiStageClass::WriteRatio()
{
   list<StageClass*>::iterator stageIter=Stages.begin();
   double prevActionChange=0.0;
   while (stageIter!=Stages.end()){
     //    cerr<<"Some stage is writing their ratio"<<endl;
     (*stageIter)->WriteRatio();
     stageIter++;
   }  
   MoveClass::WriteRatio();
   CenterOfMassVar.Write(cm2);
   cm2=0;
}


void MultiStageClass::MakeMove()
{
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  //  NewMoveProb=1.0;
  //  OldMoveProb=1.0;
  //  cerr<<"In "<<endl;
  
  struct timeval start, end;
  struct timezone tz;



  while (stageIter!=Stages.end() && toAccept){
    gettimeofday(&start, &tz);
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
  gettimeofday(&end, &tz);
  TimeSpent2 += (double)(end.tv_sec-start.tv_sec) +
    1.0e-6*(double)(end.tv_usec-start.tv_usec);
  

//     if (toAccept){
//       NewMoveProb*=(*stageIter)->NewSample*(*stageIter)->AcceptProb;
//       OldMoveProb*=(*stageIter)->OldSample*(*stageIter)->OldAcceptProb;
//     }
//     else {
//       NewMoveProb*=(*stageIter)->NewSample*(1.0-((*stageIter)->AcceptProb));
//       OldMoveProb*=(*stageIter)->OldSample*(1.0-((*stageIter)->OldAcceptProb));
//     }
//     cerr<<(*stageIter)->OldSample<<" "<<(*stageIter)->OldAcceptProb<<endl;

    stageIter++;
  }

  if (toAccept){
    Accept();
  }
  else {
    Reject();
  }
  //MoveClass::MakeMove();
  //  cerr<<NewMoveProb<<" "<<OldMoveProb<<endl;
}


void MultiStageClass::Accept()
{
  //  cerr<<"going to accept"<<endl;
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();
       stageIter!=Stages.end();stageIter++){
    (*stageIter)->Accept();
  }  
  NumAccepted++;
  cm2=cm2+PathData.Path.cm2;
  
}

void MultiStageClass::Reject()
{
  //  cerr<<"Going to reject"<<endl;
  PathData.RejectMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator 
	 stageIter=Stages.begin();stageIter!=Stages.end();stageIter++){
    (*stageIter)->Reject();
  }
}

