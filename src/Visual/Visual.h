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
  vector<Vec3> Path;
  vector<Vec3> Color;
  bool Closed;
};

class BoxClass
{
private:
  Vec3 Box, BoxInv;
  inline void PutInBox(Vec3 &r);
public:
  inline void Set (Vec3 box) 
  { 
    Box = box; 
    BoxInv[0]=1.0/Box[0]; BoxInv[1]=1.0/Box[1]; BoxInv[2]=1.0/Box[2];
  }
  inline void Set (double lx, double ly, double lz)
  {
    Box[0] = lx; Box[1] = ly; Box[2] = lz;
    BoxInv[0]=1.0/lx; BoxInv[1]=1.0/ly; BoxInv[2]=1.0/lz;
  }
  inline double operator[](int i) const
  { return Box[i]; }
  inline double& operator[](int i)
  { return Box[i]; }
  inline operator Vec3() const
  { return Box; }

  bool BreakSegment (Vec3 &r1, Vec3 &r2, Vec3 &wall1, Vec3 &wall2);

  void PutPathsInBox (vector<OnePath*>& inList,
		      vector<OnePath*>& outList);
};

inline void BoxClass::PutInBox (Vec3 &r)
{
  for (int i=0; i<3; i++) {
    double n = -floor(r[i]*BoxInv[i]+0.5);
    r[i] += n*Box[i];
  }
}

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
