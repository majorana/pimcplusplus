#ifndef MOVECHOICE_H
#define MOVECHOICE_H

/// An abstract class that lets you choose which move to make in an
/// arbitrary way 
class MoveChoiceClass
{
  /// Function to return the integer index of a move
  virtual int ChooseMove()=0;
};



/// inherits from MoveChoiceClass. The moves are chosen
/// stochastically. Broken class. Don't assume it does this.
class StochasticChoiceClass:public MoveChoiceClass
{
  /// Function to return the integer index of a move (stochastically)
  int ChooseMove();
};

///inherits from MoveChoiceClass. The moves are chosen by iterating
///over all possiblemoves. Broken class. Don't assume it does this.
class LoopChoiceClass:public MoveChoiceClass
{
 private:
  int BisectionNum = 0;
  inline int ChooseMove();
  int NumTimeSlices;


};

/// Note: THIS DOES NOT WORK CORRECTLY EITHER
inline  int  LoopChoiceClass::ChooseMove()
{
  BisectionNum++; 
  if(BisectionNum==NumTimeSlices)
    return 0;
  else
    return 1; 

}
/// Note: THIS DOES NOT WORK CORRECTLY!!!!
inline  int StochasticChoiceClass::ChooseMove()
{
  return 0;///so far first move in the list
}

#endif
