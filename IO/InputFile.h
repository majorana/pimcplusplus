#ifndef INPUTFILE_H
#define INPUTFILE_H
#include "Blitz.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stream.h>
#include <string>

inline void Abort (char *str)
{
  cerr << str << "\nAborting.\n";
  exit(1);
}

class InputBuffer
{
private:
  int pos;
  Array<char,1> buffer;
  int size;
  int StartString (char *str);

public:

  char operator()(int i) const
  {
    return (buffer(i));
  }

  void Resize(int newsize)
  {
    size = newsize;
    buffer.resize(newsize);
  }
  
  void Rewind()
  {
    pos = 0;
  }

  int Read (char *FileName);
  void Write (FILE *fout);
  inline int Size()
  {
    return (size);
  }
  int FindBlock (char StartChar, char EndChar, InputBuffer &BlockBuffer);
  int FindQuoteBlock (InputBuffer &BlockBuffer);
  int FindName (char *Name);
  int FindSection (char *SecName, InputBuffer &SectionBuff, bool rewind=true);
  int FindVarBuf (char *VarName, InputBuffer &SectionBuff);

  int ReadDouble (double &num);
  int ReadInt (int &num);
  int ReadBool (bool &IsTrue);
  int ReadQuoteString(string &str);
  int ReadVector(Array<double,1> &vec);
  int ReadVector(Array<int,1> &vec);
  int ReadVector(Array<string,1> &vec);
  int ReadString(char str[], int max);
  int ReadVar (char *VarName, double &num, bool rewind=true);
  int ReadVar (char *VarName, int &num, bool rewind=true);
  int ReadVar (char *VarName, Array<scalar,1> &vec, bool rewind=true);
  int ReadVar (char *VarName, Array<int,1> &vec, bool rewind=true);
  int ReadVar (char *VArName, Array<string,1> &vec, bool rewind=true);
  int ReadVar (char *VarName, char *str, int maxlength, bool rewind=true);
  int ReadVar (char *VarName, bool &IsTrue, bool rewind=true);
  
  InputBuffer()
  {
    size = 0;
    pos = 0;
  }
};



#endif

