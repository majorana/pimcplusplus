#ifndef MOVECHOICE_H
#define MOVECHOICE_H

///abstract class
class MoveChoiceClass
{
  virtual int ChooseMove()=0;
};



///inherits from MoveChoiceClass
class StochasticChoiceClass:public MoveChoiceClass
{
  int ChooseMove();
};

///inherits from MoveChoiceClass
class LoopChoiceClass:public MoveChoiceClass
{
 private:
  int BisectionNum = 0;
  inline int ChooseMove();
  int NumTimeSlices;


};

inline  int  LoopChoiceClass::ChooseMove()
{
  BisectionNum++; 
  if(BisectionNum==NumTimeSlices)
    return 0;
  else
    return 1; 

}

inline  int StochasticChoiceClass::ChooseMove()
{
  return 0;//so far first move in the list
}

#endif
