#ifndef VISUAL_H
#define VISUAL_H

#include "PathVis.h"
#include "../Common/IO/InputOutput.h"
#include <gtkmm/adjustment.h>


class SpeciesClass
{
public:
  double lambda;
  string Name;
  int NumParticles;
  int FirstParticle, LastParticle;
  SpeciesClass () : FirstParticle(0)
  { 

  }
};

typedef enum {LINES, TUBES} PathTypeType;

class VisualClass : public Gtk::Window
{
protected:
  // signal handlers:
  void on_button_quit_clicked();

  Array<double,4> Paths;
  Array<SpeciesClass,1> Species;
  Vec3 Box;

  // member widgets:
  Gtk::VBox m_VBox;
  Gtk::Button m_ButtonQuit;
  Gtk::HScale FrameScale;
  Gtk::Adjustment FrameAdjust;
  Gtk::Toolbar Tools;
  Gtk::RadioToolButton LinesButton, TubesButton, SmoothButton, StraightButton;
  void FrameChanged();
  PathTypeType PathType; 
  void LineToggle();
  Gtk::Image TubesImage, LinesImage, StraightImage, SmoothImage;
  Gtk::SeparatorToolItem ToolSep;

  Glib::RefPtr<Gtk::ActionGroup> Actions;
  Glib::RefPtr<Gtk::UIManager> Manager;

  void OnOpen();
  void OnExport();
  //  bool on_delete_event();
  void Quit();
  void ResetView();
  void PutInBox();

  Gtk::FileChooserDialog FileChooser;
public:
  PathVisClass PathVis;

  void Read(string fileName);
  void MakeFrame (int frame);
  
  VisualClass();
  virtual ~VisualClass();

};

#endif
