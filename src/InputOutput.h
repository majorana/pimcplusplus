#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <string>
#include <list>

using namespace std;


class InputSectionClass
{
private:
  list<InputSectionClass*> MyChildren;
  InputSection *MyParent;
public:
  virtual void Open (char *FileName) = 0;
  void Open (string FileName) { Open (FileName.c_str());}

  InputSectionClass& FindSection (string Name, bool rewind=true);
  /// Atomics
  virtual bool ReadVar (string name, int var)    = 0;
  virtual bool ReadVar (string name, double var) = 0;
  virtual bool ReadVar (string name, bool var)   = 0;
  virtual bool ReadVar (string name, string var) = 0;
  /// 1D Arrays
  virtual bool ReadVar (string name, Array<int,1> var)    = 0;
  virtual bool ReadVar (string name, Array<double,1> var) = 0;
  virtual bool ReadVar (string name, Array<bool,1> var)   = 0;
  virtual bool ReadVar (string name, Array<bool,1> var)   = 0;
  /// 1D Arrays
  virtual bool ReadVar (string name, Array<int,2> var)    = 0;
  virtual bool ReadVar (string name, Array<double,2> var) = 0;
  virtual bool ReadVar (string name, Array<bool,2> var)   = 0;
  virtual bool ReadVar (string name, Array<bool,2> var)   = 0;

  template<class T> bool ReadVar (string name, Array<T,1> var);

};


class OutputSectionClass
{
public:
  


};



class Object
{
public:
  virtual string Name = 0;
  virtual Read (InputSectionClass &section) = 0;
  virtual CheckPoint(OutputSectionClass &outSection);
};



#endif
