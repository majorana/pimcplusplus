#include "PathVis.h"

namespace Trackball {
  extern "C" {
    #include "trackball.h"
  }
}

ViewClass::ViewClass (PathVisClass &pathVis) :
  PathVis(pathVis), Button1Pressed(false), MinScale(0.2), MaxScale(5.0),
  Scale(1.0), Distance(3.0), UsePerspective(false)
{
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      RotMat[i][j] = (i==j) ? 1.0 : 0.0;
  Quaternion[0] = 0.0;  Quaternion[1] = 0.0;
  Quaternion[2] = 0.0;  Quaternion[3] = 1.0;
  PathVis.signal_button_press_event().connect
    (sigc::mem_fun(*this, &ViewClass::OnButtonPress));
  PathVis.signal_button_release_event().connect
    (sigc::mem_fun(*this, &ViewClass::OnButtonRelease));
  PathVis.signal_motion_notify_event().connect
    (sigc::mem_fun(*this, &ViewClass::OnMotion));

  PathVis.add_events(Gdk::BUTTON1_MOTION_MASK    |
		     Gdk::BUTTON2_MOTION_MASK    |
		     Gdk::BUTTON_PRESS_MASK      |
		     Gdk::BUTTON_RELEASE_MASK    |
		     Gdk::VISIBILITY_NOTIFY_MASK);
}


bool 
ViewClass::OnButtonPress (GdkEventButton *event)
{
  if (event->button == 1) {
    StartX = event->x;
    StartY = event->y;
    Button1Pressed=true;
    Button2Pressed=false;
  }
  else if (event->button == 2) {
    StartX = event->x;
    StartY = event->y;
    Button1Pressed=false;
    Button2Pressed=true;
    OldScale = Scale;
  }
  return false;
}


bool 
ViewClass::OnButtonRelease (GdkEventButton *event)
{
  if (event->button == 1)
    Button1Pressed = false;
  if (event->button == 2)
    Button2Pressed = false;
  return false;
}

bool 
ViewClass::OnMotion (GdkEventMotion *event)
{
  if (Button1Pressed) {
    double w = PathVis.get_width();
    double h = PathVis.get_height();
    double x = event->x;
    double y = event->y;
    // scale 
    double uStart = (2.0*StartX - w) / w;
    double vStart = (h - 2.0*StartY) / h;
    double uEnd   = (2.0*x      - w) / w;
    double vEnd   = (h - 2.0*y     ) / h;
    
    double dQuat[4];
    
    Trackball::trackball (dQuat, uStart, vStart, uEnd, vEnd);
    Trackball::add_quats (dQuat, Quaternion, Quaternion);
    Trackball::build_rotmatrix (RotMat, Quaternion);
    PathVis.Invalidate();
    StartX = x;
    StartY = y;
  }
  else if (Button2Pressed) {
    double w = PathVis.get_width();
    double delta = (event->x - StartX) / w;
    if (delta < 0.0)
      Scale = max (MinScale, OldScale + delta);
    else
      Scale = min (MaxScale, OldScale + delta);
    PathVis.Invalidate();
  }
  return false;
}


void 
ViewClass::GLtransform()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //  gluPerspective(40.0, 1.0, 1.0, 10.0);
  if (UsePerspective)
    gluPerspective(40.0, 1.0, 1.0, 8.0*Distance/Scale);
  else
    glOrtho(-0.6*Distance/Scale, 0.6*Distance/Scale, 
	    -0.6*Distance/Scale, 0.6*Distance/Scale, 1.0, 8.0*Distance/Scale);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(0.0, 0.0, Distance/Scale,
            0.0, 0.0, 0.0,
            0.0, 1.0, 0.0);
  
  glTranslatef(0.0, 0.0, -Distance/Scale);
  //  glScaled(Scale, Scale, Scale);
  glMultMatrixd(&RotMat[0][0]);
}


void 
ViewClass::POVtransform (FILE *fout)
{
  fprintf (fout, "camera {\n");
  fprintf (fout, "  location <%14.10f, %14.10f %14.10f>\n",
	   0.0, 0.0, 2.0*Distance/Scale);
  double angle = 51.661;
  angle = 40.0;
  fprintf (fout, "  angle %1.5f\n", angle);
  fprintf (fout, "  right <-1.0,0,0>\n");
  fprintf (fout, "  look_at <%14.10f %14.10f %14.10f>\n}\n\n",
	   0.0, 0.0, 0.0);
}


void 
ViewClass::SetDistance (double dist)
{
  Distance = dist;
}


void 
ViewClass::Reset()
{
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      RotMat[i][j] = (i==j) ? 1.0 : 0.0;
  Quaternion[0] = 0.0;  Quaternion[1] = 0.0;
  Quaternion[2] = 0.0;  Quaternion[3] = 1.0;
  Scale = 1.0;
}


string
ViewClass::RotationString()
{
  char rotString[500];
  
  snprintf(rotString,500,
	   "  matrix <%8.5f, %8.5f, %8.5f,\n          %8.5f, %8.5f, %8.5f,\n          %8.5f, %8.5f, %8.5f,\n          %8.5f, %8.5f, %8.5f >\n",
	    RotMat[0][0], RotMat[0][1], RotMat[0][2],
	    RotMat[1][0], RotMat[1][1], RotMat[1][2],
	    RotMat[2][0], RotMat[2][1], RotMat[2][2],
	    0.0         , 0.0         , 0.0         );
  string str = rotString;
  return str;
}
