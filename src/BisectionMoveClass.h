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




void getNextParticles(Array<int,1> &theParticles)
{

  for (int counter=0;counter<theParticles.size()){    
    bool particleOk=false;
    while (!particleOk){
      theParticle(counter)=floor(sprng()*NumOfParticlesToMove);
      particleOk=true;
      for (int nter2=0;counter2<counter;counter2++){
	if (theParticle(counter)==theParticle(counter2)){
	  particleOk=false;
	}
      }
    }
  }

}

void BisectionMoveClass::makeMove()
{

  Array<int,1> theParticles;
  theParticles.resize(NumOfParticlesToMove);
    
  int EndTimeSlice=1<<NumLevels+StartSlice;
  skip=1<<NumLevels;
  for (int levelCounter=NumLevels;levelCounter>0;NumLevels++){
    for (int Slice=StartTimeSlice+skip;Slice<EndTimeSlice;Slice=Slice+skip){
      getNextParticles(theParticles);
      (*PathData).


      singleSliceSample(theParticles,Slice);
    }

    skip=skip>>1
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
