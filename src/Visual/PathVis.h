#include "../Common/IO/InputOutput.h"
#include <iostream>
#include <cstdlib>

#include <gtkmm.h>

#include <gtkglmm.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "../Common/Blitz.h"
#include "GLObject.h"
#include "PathObject.h"
#include "BoxObject.h"
#include "SphereObject.h"


class PathVisClass;

class ViewClass : public sigc::trackable
{
private:
  friend class PathVisClass;
  double StartX, StartY;
  bool Button1Pressed, Button2Pressed; 
  double MinScale, MaxScale;
  PathVisClass &PathVis;
  double Distance;
  bool UsePerspective;
public:
  double Scale, OldScale;
  double Quaternion[4];
  double RotMat[4][4];

  bool OnButtonPress   (GdkEventButton* event);
  bool OnButtonRelease (GdkEventButton* event);
  bool OnMotion        (GdkEventMotion* event);
  void SetDistance (double dist);

  inline bool SetPerspective (bool usePersp) 
  { UsePerspective = usePersp; }
  void GLtransform();
  void Reset();

  ViewClass (PathVisClass &pathVis);
};

class PathVisClass : public Gtk::DrawingArea,
		     public Gtk::GL::Widget<PathVisClass>
{
  friend class ViewClass;
protected:
  Glib::RefPtr<Gdk::GL::Window> GLwindow;

  virtual void on_realize();
  virtual bool on_configure_event(GdkEventConfigure* event);
  virtual bool on_expose_event(GdkEventExpose* event);
  int NumLists;
public:
  ViewClass View;
  vector<GLObject *> Objects;
  void AddBox  (double xSize, double ySize, double zSize);
  void AddPath (Array<Vec3,1> &path, bool closed=true);
  void GLRender();
  void Invalidate();

  // Constructor
  PathVisClass();
  // Destructor
  virtual ~PathVisClass();
};

