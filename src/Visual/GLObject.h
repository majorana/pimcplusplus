#ifndef GL_OBJECT_H
#define GL_OBJECT_H
#include <stdio.h>

class GLObject
{
protected:
  int ListNum;
  bool Created;
  bool Visible;
public:
  virtual void Draw();
  virtual void DrawPOV(FILE *fout) = 0;
  void Start();
  void End();
  void SetVisible (bool visible);
  bool GetVisible ();
  GLObject();
  ~GLObject();
};

#endif
