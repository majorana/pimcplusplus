#ifndef INPUT_OUTPUT_ASCII_H
#define INPUT_OUTPUT_ASCII_H


#include "InputOutputBase.h"



class InputSectionASCIIClass : public InputSectionClass
{

  void ReadWithoutComments(string fileName, Array<char,1> &buffer);

 public:
  void PrintTree(InputSectionClass *sec);
  void PrintTree() { }
  void PrintTree(int index) { }
  bool OpenFile (string fileName, InputSectionClass *parent);
  bool OpenFile (string filename, string parentName, 
		 InputSectionClass *parent) { }
  void Close() { };
  void CloseFile();

  
  
 
};


class VarASCIIClass : public VarClass
{  
public:
  Array<char,1> Value;
  bool ReadInto (double &val)
  {
    int l = Value.size();
    char str[l+1];
    for (int i=0; i<l; i++)
      str[i] = Value(i);
    str[l] ='\0';
    val = atof(str);
    return true;
  }
  bool ReadInto (int &val){
    int l = Value.size();
    char str[l+1];
    for (int i=0; i<l; i++)
      str[i] = Value(i);
    str[l] ='\0';
    cerr<<"Value is "<<Value<<endl;
    val = atoi(str);
    cerr<<"val is "<<val<<endl;
    return true;
  }
  bool ReadInto (string &val){
    val="";
    for (int counter=1;counter<Value.size()-1;counter++){
      val=val+Value(counter);
    }
    return true;
  }
  bool ReadInto (Array<double,1> &v) { }
  bool ReadInto (Array<double,2> &v) { }
  bool ReadInto (Array<double,3> &v) { }
  bool ReadInto (Array<int,1> &v) { }
  bool ReadInto (Array<int,2> &v) { }
  bool ReadInto (Array<int,3> &v) { }
  bool ReadInto (Array<string,1> &v) { }
  bool ReadInto (Array<string,2> &v) { }
  bool ReadInto (Array<string,3> &v) { }

};

#endif
