#ifndef BISECTIONMOVE_CLASS_H
#define BISECTIONMOVE_CLASS__H

#include "MoveClass.h";


class BisectionMoveClass::MoveClass
{
private:

  int NumOfParticlesToMove;
  int ParticleTypeToMove;
  int StartTimeSlice;
  int NumLevels;
  
  void makeMove();
 
};




void BisectionMoveClass::makeMove()
{
  bool accept=true;
  double oldLogSampleProb;
  double newLogSampleProb;
  Array<int,1> theParticles;
  double logSampleProb;
  int EndTimeSlice=1<<NumLevels+StartSlice;

  for (int levelCounter=NumLevels;levelCounter>0;NumLevels++){
    ChooseParticles();   
    setMode(OLD);

    oldAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,level);
    setMode(NEW);
    (*PathData).TotalAction.SampleParticles(ActiveParticles,StartSlice,EndSlice,levelCounter,newLogSampleProb,oldLogSampleProb);
    newAction = (*PathData).TotalAction.calcTotalAction(ActiveParticles,level);
    
    
    
  }
  

}


       
  


  int startLoc;
  int endLoc;
  int timeSliceToMove;
  thePathToMove=choosePathToMove(pathData);
  chooseSegmentToMove(pathData,&startLoc,&endLoc);
  segmentSpan=endLoc-startLoc; 


 
  for (int denominator=2;denominator=denominator*2;denominator<=segmentSpan){  //this is a loop over the level that the bisection move is on
    for (int numerator=1;numerator<denominator;denominator=denominator+2){//this is a loop over all the particles in that specific level
      timeSliceToMove=numerator/denominator*segmentSpan
	thePath.storeState(timeSliceToMove,thePathToMove,pathData);//I think these are sent as pointers of objects but I'm not positive about that.
      singleLinkMove(timeSliceToMove,thePathToMove,pathData);


      //here we want to add to some list that we
      //have made a change in this specific link so
      //somebody knows how to deal with it.
    }
    //here we should go about checking the action at this
    //level and deciding if we want to accept or reject
    //this level.
  }
			  














#endif
