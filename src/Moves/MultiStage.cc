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

#include "MultiStage.h"


void MultiStageClass::Read(IOSectionClass& in)
{
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
}


void MultiStageClass::MakeMove()
{
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  //  NewMoveProb=1.0;
  //  OldMoveProb=1.0;
  cerr<<"In "<<endl;
  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);

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
  //  cerr<<"Out "<<endl;
  if (toAccept)
    Accept();
  else 
    Reject();
  //MoveClass::MakeMove();
  //  cerr<<NewMoveProb<<" "<<OldMoveProb<<endl;
}


void MultiStageClass::Accept()
{
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();
       stageIter!=Stages.end();stageIter++){
    (*stageIter)->Accept();
  }  
  NumAccepted++;
}

void MultiStageClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator 
	 stageIter=Stages.begin();stageIter!=Stages.end();stageIter++){
    (*stageIter)->Reject();
  }
}

