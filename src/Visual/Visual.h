#ifndef VISUAL_H
#define VISUAL_H

#include "PathVis.h"
#include "../Common/IO/InputOutput.h"
#include <gtkmm/adjustment.h>
#include "OnePath.h"
#include "BoxClass.h"
#include "SmoothClass.h"

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
  //////////
  // Data //
  //////////
  Array<double,4> PathArray;
  Array<int,2> PermArray;
  Array<SpeciesClass,1> Species;
  vector<OnePath*> Paths;
  BoxClass Box;
  SmoothClass Smoother;
  void MakePaths(int frame);

  // member widgets:
  Gtk::VBox m_VBox;
  Gtk::Button m_ButtonQuit;
  Gtk::HScale FrameScale;
  Gtk::Adjustment FrameAdjust;
  Gtk::Toolbar Tools;
  Gtk::HBox ToolBox;
  Gtk::RadioToolButton LinesButton, TubesButton, 
    SmoothButton, StraightButton,
    WrapButton, NoWrapButton;
  void FrameChanged();
  PathTypeType PathType; 
  void LineToggle(), WrapToggle(), SmoothToggle();
  Gtk::Image TubesImage, LinesImage, StraightImage, SmoothImage,
    NoWrapImage, WrapImage;
  Gtk::SeparatorToolItem ToolSep1, ToolSep2;

  // Detail control
  Gtk::Frame DetailFrame;
  Gtk::HScale DetailScale;
  Gtk::Adjustment DetailAdjust;


  Glib::RefPtr<Gtk::ActionGroup> Actions;
  Glib::RefPtr<Gtk::UIManager> Manager;

  void OnOpen();
  void OnExport();
  //  bool on_delete_event();
  void Quit();
  void ResetView();
  void PutInBox();
  void OnDetailChange();

  Gtk::FileChooserDialog FileChooser;

  bool Wrap, Smooth;
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
