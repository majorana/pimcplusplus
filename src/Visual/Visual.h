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

class OnePath
{
public:
  Array<Vec3,1> Path;
  Array<Vec3,1> Color;
  bool Closed;
};

class VisualClass : public Gtk::Window
{
protected:
  //////////
  // Data //
  //////////
  Array<double,4> PathArray;
  Array<int,2> PermArray;
  Array<SpeciesClass,1> Species;
  vector<OnePath*> Paths;
  Vec3 Box;
  void MakePaths(int frame);

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

  inline SpeciesClass& PtclSpecies (int ptcl);
  void Read(string fileName);
  void MakeFrame (int frame);
  
  VisualClass();
  virtual ~VisualClass();

};


inline SpeciesClass& VisualClass::PtclSpecies (int ptcl)
{
  int si = 0;
  while (si < Species.size() && (Species(si).FirstParticle > ptcl))
    si++;
  return Species(si);
}

#endif
