#ifndef GL_OBJECT_H
#define GL_OBJECT_H

class GLObject
{
protected:
  int ListNum;
  bool Created;
  bool Visible;
public:
  virtual void Draw();
  void Start();
  void End();
  void SetVisible (bool visible);
  bool GetVisible ();
  GLObject();
  ~GLObject();
};

#endif
