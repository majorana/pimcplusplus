#include "WrapClass.h"

void LoopClass::Read(IOSectionClass &IO)
{
  int steps;
  assert(IO.ReadVar("steps",steps));
  Read(IO, steps);
}

void LoopClass::Read(IOSectionClass &IO, int steps)
{
  Steps = steps;
  int numSections=IO.CountSections();
  //  Array<EventClass*,1> EventList(numSections);
  Events.resize(numSections);

  for (int counter=0;counter<numSections;counter++){
    IO.OpenSection(counter);
    if (IO.GetName()=="Move"){
      cerr<<"Processing move"<<endl;
      string moveName="";
      Events(counter)=new MoveWrap();
      ((MoveWrap*)(Events(counter)))->Move=NULL;
      assert(IO.ReadVar("name",moveName));
      for (int moveCounter=0;moveCounter<(*MoveList).size();moveCounter++){
	if (moveName==(*MoveList)(moveCounter)->Name){
	  ((MoveWrap*)(Events(counter)))->Move=(*MoveList)(moveCounter);
	}
      }
      if (((MoveWrap*)(Events(counter)))->Move==NULL){
	cerr<<"We didn't find the name of the move"<<endl;
	cerr<<"It was called "<<moveName<<endl;
	assert(1==2);
      }
      cerr<<"End move"<<endl;
    }
    else if (IO.GetName()=="Observe"){
      cerr<<"Processing observable"<<endl;
      string observeName="";
      Events(counter)=new ObserveWrap();
      cerr<<"Test 1"<<endl;
      ((ObserveWrap*)(Events(counter)))->Observe=NULL;
      assert(IO.ReadVar("name",observeName));
      cerr<<"Test 2"<<endl;
      for (int obsCount=0;obsCount<(*ObserveList).size();obsCount++){
	cerr<<"Counter number is "<<obsCount<<endl;
	if (observeName==(*ObserveList)(obsCount)->Name){
	  ((ObserveWrap*)(Events(counter)))->Observe=(*ObserveList)(obsCount);
	}
      }
      cerr<<"Test 3"<<endl;
      if (((ObserveWrap*)(Events(counter)))->Observe==NULL){
	cerr<<"We didn't find the name of the observable"<<endl;
	cerr<<"It was called "<<observeName<<endl;
	for (counter=0;counter<(*ObserveList).size();counter++){
	  cerr<<"The Name is "<<(*ObserveList)(counter)->Name<<endl;
	}
	cerr<<"We did the loop"<<endl;
	assert(1==2);
      }
      cerr<<"End Observable"<<endl;
    }
    else if (IO.GetName()=="Loop"){
      cerr<<"Starting Loop"<<endl;
      cerr<<"Size is "<<(*ObserveList).size()<<endl;
      LoopClass *tempLoop=new LoopClass(MoveList,ObserveList);
      tempLoop->Read(IO);
      Events(counter)=tempLoop;
      cerr<<"Ending Loop"<<endl;
    }
    else {
      cerr<<"Don't recognize this construct"<<endl;
    }
    IO.CloseSection();
  }
}

  
