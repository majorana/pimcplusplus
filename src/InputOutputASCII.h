#ifndef INPUT_OUTPUT_ASCII_H
#define INPUT_OUTPUT_ASCII_H


#include "InputOutput.h"


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
};

#endif
