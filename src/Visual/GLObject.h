#ifndef GL_OBJECT_H
#define GL_OBJECT_H
#include <stdio.h>
#include <string>

using namespace std;

class GLObject
{
protected:
  int ListNum;
  bool Created;
  bool Visible;
  bool Transparent;
public:
  virtual void Draw();
  virtual void DrawPOV(FILE *fout, string rotString) = 0;
  void Start();
  void End();
  void SetVisible (bool visible);
  bool GetVisible ();
  void SetTransparent (bool trans);
  bool GetTransparent ();
  GLObject();
  ~GLObject();
};

#endif
