#ifndef VISUAL_H
#define VISUAL_H

#include "PathVis.h"
#include "../Common/IO/InputOutput.h"
#include <gtkmm/adjustment.h>
#include "OnePath.h"
#include "BoxClass.h"
#include "SmoothClass.h"
#include "Export.h"

/// This species class stores info about the species for the path
/// visualization.  Not the same as the pimc++ version.
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


/// VisualClass is the main GUI interaction for the path visualization
/// code.  It allows the user to rotate around and zoom in and out of
/// a box containing paths and classical particles dumped from
/// pimc++.  In includes features such as a frame slider, line or tube
/// rendering of paths, and optional fourier smoothing with detail
/// control.  
class VisualClass : public Gtk::Window
{
protected:
  friend class ExportClass;
  //////////
  // Data //
  //////////
  /// This stores the raw paths from the pimc++ output file.  The 4
  /// dimension are (frame, ptcl, slice, dim)
  Array<double,4> PathArray;
  /// This stores the global permutations for the paths.
  /// (frame, ptcl).  The permutation acts between the last and first
  /// slice.  
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
  ExportClass Export;
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
