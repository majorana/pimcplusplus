#ifndef FILE_EXPAND_H
#define FILE_EXPAND_H

#ifndef MAC
#include <wordexp.h>
#endif
#include <stdlib.h>

inline string ExpandFileName(string fname)
{
  string outName;
#ifdef MAC
  if (fname[0] == '~') {
    outName = getenv ("HOME"); 
    fname = fname.erase(0,1);
  }
  outName.append(fname);
#else
  wordexp_t words;
  wordexp (fname.c_str(), &words, 0);
  outName = words.we_wordv[0];
  wordfree(&words);
#endif
  return outName;
}

#endif
