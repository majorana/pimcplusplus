#ifndef INPUT_OUTPUT_ASCII_H
#define INPUT_OUTPUT_ASCII_H
#include "InputOutputBase.h"
#include <iostream>
#include <stack>
#include <string>
#include <list>



/// This class holds an ASCII token, which is just a string and the
/// line number in which it appeared in the file.
class TokenClass
{
public:
  string Str;
  int LineNumber;
};


/// This is the ASCII specialization of InputTreeClass for ASCII text
/// files.  It's syntax is as follows:
/// Section (SectionName)
/// {
///   double x = 3;
///   Array<int,1> y(3) = [1, 2, 3];
///   Array<int,3> z(2,2,1) = [ 1, 2, 
///                             3, 4 ];
///   Section (cabbage, "cabbagesection.txt")
/// }
class InputTreeASCIIClass : public InputTreeClass
{

  void ReadWithoutComments(string fileName, Array<char,1> &buffer);
  bool ReadSection (InputTreeClass *parent, string name,
		    list<TokenClass>::iterator &iter,
		    list<TokenClass> &tokenList,
		    bool wantEndBrace);

 public:
  void PrintTree(int level);
  void PrintTree();
  bool OpenFile (string filename, string parentName, 
		 InputTreeClass *parent);
  void Close() { };
  void CloseFile();

  
  
 
};


class VarASCIIClass : public VarClass
{  
public:
  void *Value;

  bool ReadInto (double &val);
  bool ReadInto (int &val);
  bool ReadInto (string &val);
  bool ReadInto (bool &val);
  bool ReadInto (Array<double,1> &v);
  bool ReadInto (Array<double,2> &v);
  bool ReadInto (Array<double,3> &v);
  bool ReadInto (Array<int,1> &v);
  bool ReadInto (Array<int,2> &v);
  bool ReadInto (Array<int,3> &v);
  bool ReadInto (Array<string,1> &v);
  bool ReadInto (Array<string,2> &v);
  bool ReadInto (Array<string,3> &v);
  bool ReadInto (Array<bool,1> &v);
  bool ReadInto (Array<bool,2> &v);
  bool ReadInto (Array<bool,3> &v);


};

#endif
