#include "Visual.h"

void VisualClass::MakeFrame(int frame)
{
  int numSpecies = Species.size();

  for (int i=0; i<PathVis.Objects.size(); i++)
    delete PathVis.Objects[i];
  
  PathVis.Objects.resize(0);

  MakePaths(frame);

  int numPtcls = PathArray.extent(1);
  int numSlices = PathArray.extent(2);


//   Array<Vec3,1> onePath(numSlices);
//   for (int si=0; si<numSpecies; si++) {
//     for (int ptcl=Species(si).FirstParticle; 
// 	 ptcl<=Species(si).LastParticle; ptcl++) {
//       if (Species(si).lambda != 0.0) {
// 	for (int slice=0; slice<numSlices; slice++) {
// 	  onePath(slice)[0] = PathArray(frame,ptcl,slice,0);
// 	  onePath(slice)[1] = PathArray(frame,ptcl,slice,1);
// 	  onePath(slice)[2] = PathArray(frame,ptcl,slice,2);
// 	}
// 	PathObject* pathObj = new PathObject;
// 	if (si == 0)
// 	  pathObj->SetColor (0.3, 0.3, 1.0);
// 	else
// 	  pathObj->SetColor (0.0, 1.0, 0.0);
// 	if (PathType == TUBES)
// 	  pathObj->TubesSet (onePath);
// 	else
// 	  pathObj->LinesSet (onePath);
// 	PathVis.Objects.push_back(pathObj);
//       }
//       else {
// 	Vec3 pos;
// 	pos[0] = PathArray(frame, ptcl, 0, 0);
// 	pos[1] = PathArray(frame, ptcl, 0, 1);
// 	pos[2] = PathArray(frame, ptcl, 0, 2);

// 	SphereObject* sphere = new SphereObject;
// 	sphere->SetPos (pos);
// 	sphere->SetColor (Vec3(1.0, 0.0, 1.0));
// 	PathVis.Objects.push_back(sphere);
//       }
//     }
//   }
  cerr << "Paths.size() = " << Paths.size() << endl;
  for (int li=0; li<Paths.size(); li++) {
    PathObject* pathObj = new PathObject;
    pathObj->SetColor (0.0, 0.0, 1.0);
    pathObj->SetRadius (min(min(Box[0], Box[1]), Box[2])*0.005);
    if (PathType == TUBES)
      pathObj->TubesSet (Paths[li]->Path);
    else
      pathObj->LinesSet (Paths[li]->Path);
    PathVis.Objects.push_back(pathObj);
  }
  for (int si=0; si<numSpecies; si++) 
    if (Species(si).lambda == 0)
      for (int ptcl=Species(si).FirstParticle; ptcl<=Species(si).LastParticle;
	   ptcl++) {
	SphereObject* sphere = new SphereObject;
	dVec pos;
	pos[0] = PathArray(frame, ptcl, 0, 0);
 	pos[1] = PathArray(frame, ptcl, 0, 1);
 	pos[2] = PathArray(frame, ptcl, 0, 2);
	sphere->SetPos (pos);
	sphere->SetColor (Vec3(1.0, 0.0, 1.0));
	PathVis.Objects.push_back(sphere);
      }
  
  BoxObject *boxObject = new BoxObject;
  boxObject->SetColor (0.5, 0.5, 1.0);
  boxObject->Set (Box[0], Box[1], Box[2]);
  PathVis.Objects.push_back(boxObject);
}

void VisualClass::Read(string fileName)
{
  IOSectionClass in;
  assert(in.OpenFile (fileName));

  Array<double,1> box;
  assert (in.OpenSection("System"));
  assert (in.ReadVar ("Box", box));
  Box.Set (box(0), box(1), box(2));
  cerr << "Box = " << box << endl;

  double maxDim = max(max(box(0), box(1)), box(2));
  PathVis.View.SetDistance (1.2*maxDim);
  //PathVis.View.SetDistance (0.2*maxDim);

  int numSpecies = in.CountSections ("Species");  
  Species.resize(numSpecies);
  for (int i=0; i<numSpecies; i++)
    Species(i).FirstParticle = 0;
  for (int i=0; i<numSpecies; i++) {
    in.OpenSection("Species",i);
    assert (in.ReadVar("lambda", Species(i).lambda));
    assert (in.ReadVar("Name", Species(i).Name));
    assert (in.ReadVar("NumParticles", Species(i).NumParticles));
    for (int j=i+1; j<numSpecies; j++)
      Species(j).FirstParticle += Species(i).NumParticles;
    Species(i).LastParticle=Species(i).FirstParticle+Species(i).NumParticles-1;
    in.CloseSection(); // Species
  }

  for (int i=0; i<numSpecies; i++) 
    cerr << "Species:  Name = " << Species(i).Name << "    First Ptcl = " 
	 << Species(i).FirstParticle 
	 << "   Last Ptcl = " << Species(i).LastParticle << "\n";

  in.CloseSection (); // "System"

  assert(in.OpenSection("Observables"));
  assert(in.OpenSection("PathDump"));
  assert(in.ReadVar ("Path", PathArray));
  assert(in.ReadVar ("Permutation", PermArray));
  PutInBox();

  FrameAdjust.set_upper(PathArray.extent(0)-1);
  
  in.CloseSection();
  in.CloseFile();
  FrameScale.set_value(0.0);
  MakeFrame (0);
}

void VisualClass::PutInBox()
{
  for (int frame=0; frame<PathArray.extent(0); frame++) 
    for (int ptcl=0; ptcl<PathArray.extent(1); ptcl++) 
      for (int dim=0; dim<3; dim++) {
	while (PathArray(frame,ptcl,0,dim) > 0.5*Box[dim])
	  for (int slice=0; slice<PathArray.extent(2); slice++)
	    PathArray(frame,ptcl,slice,dim) -= Box[dim];
	while (PathArray(frame,ptcl,0,dim) < -0.5*Box[dim])
	  for (int slice=0; slice<PathArray.extent(2); slice++)
	    PathArray(frame,ptcl,slice,dim) += Box[dim];
      }
	  
}


#include "tubes.xpm"
#include "lines.xpm"

VisualClass::VisualClass()
  : m_VBox(false, 0), m_ButtonQuit("Quit"), 
    FrameAdjust (0.0, 0.0, 0.0),
    PathType (LINES),
    TubesImage("tubes.png"), LinesImage("lines.png"),
    StraightImage("straight.png"), SmoothImage("smooth.png"),
    FileChooser ("Choose an output file"),
    Wrap(false)
{
  // Top-level window.
  set_title("VisualClass");

  // Get automatically redrawn if any of their children changed allocation.
  set_reallocate_redraws(true);
  add(m_VBox);

  // VisualClass OpenGL scene.
  PathVis.set_size_request(600, 600);


  // VisualClass quit button.

  m_ButtonQuit.signal_clicked().connect
    (sigc::mem_fun(*this, &VisualClass::Quit));
  FrameScale.set_adjustment (FrameAdjust);
  FrameScale.signal_value_changed().connect
    (sigc::mem_fun(*this, &VisualClass::FrameChanged));
  FrameAdjust.set_step_increment(1.0);
  FrameScale.set_digits(0);

  // Setup tool bar
  Gtk::RadioButtonGroup group = LinesButton.get_group();
  TubesButton.set_group (group);
  group = StraightButton.get_group();
  SmoothButton.set_group(group);
  group = NoWrapButton.get_group();
  WrapButton.set_group(group);

  LinesButton.set_label("Lines");
  TubesButton.set_label("Tubes");
  StraightButton.set_label("Straight");
  SmoothButton.set_label("Smooth");
  NoWrapButton.set_label("No Wrap");
  WrapButton.set_label("Wrap");

  LinesButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::LineToggle));
  WrapButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::WrapToggle));

  TubesButton.set_icon_widget (TubesImage);
  LinesButton.set_icon_widget (LinesImage);
  StraightButton.set_icon_widget(StraightImage);
  SmoothButton.set_icon_widget(SmoothImage);
  Tools.append (LinesButton);
  Tools.append (TubesButton);
  Tools.append (ToolSep1);
  Tools.append (StraightButton);
  Tools.append (SmoothButton);
  Tools.append (ToolSep2);
  Tools.append (NoWrapButton);
  Tools.append (WrapButton);

  // Setup the file chooser
  FileChooser.set_select_multiple(false);
  FileChooser.set_action(Gtk::FILE_CHOOSER_ACTION_OPEN);
  FileChooser.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  FileChooser.add_button("OK", Gtk::RESPONSE_OK);


//   Gtk::Color black ("Black");
//   Glib::RefPtr<Pixmap> linesPM = 
//     Gdk::Pixmap::create_from_xpm (TubesButton.window, NULL, black, lines);


  Actions = Gtk::ActionGroup::create();
  Actions->add (Gtk::Action::create("MenuFile", "_File"));
  Actions->add (Gtk::Action::create("Open", "_Open"),
		sigc::mem_fun(*this, &VisualClass::OnOpen));
  Actions->add (Gtk::Action::create("Export", "_Export POV"),
		sigc::mem_fun(*this, &VisualClass::OnExport));
  Actions->add (Gtk::Action::create("Quit", "_Quit"),
		sigc::mem_fun(*this, &VisualClass::Quit));
  Actions->add (Gtk::Action::create("MenuView", "View"));
  Actions->add (Gtk::Action::create("Reset", "Reset"),
		sigc::mem_fun(*this, &VisualClass::ResetView));

  Manager = Gtk::UIManager::create();
  Manager->insert_action_group(Actions);
  add_accel_group (Manager->get_accel_group());

  Glib::ustring ui_info =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='MenuFile'>"
    "      <menuitem action='Open'/>"
    "      <menuitem action='Export'/>"
    "      <separator/>"
    "      <menuitem action='Quit'/>"
    "    </menu>"
    "    <menu action='MenuView'>"
    "      <menuitem action='Reset'/>"
    "    </menu>"
    "  </menubar>"
    "  <toolbar  name='ToolBar'>"
    "    <toolitem action='Open'/>"
    "    <toolitem action='Quit'/>"
    "  </toolbar>"
    "</ui>";
  
  Manager->add_ui_from_string (ui_info);
  m_VBox.pack_start (*Manager->get_widget("/MenuBar"), Gtk::PACK_SHRINK, 0);
  m_VBox.pack_start(Tools, Gtk::PACK_SHRINK, 0);
  m_VBox.pack_start(PathVis);
  m_VBox.pack_start(FrameScale, Gtk::PACK_SHRINK,0);
  m_VBox.pack_start(m_ButtonQuit, Gtk::PACK_SHRINK, 0);

  // Show window.
  show_all();
}

VisualClass::~VisualClass()
{}
  
void VisualClass::OnOpen()
{
  int result = FileChooser.run();
  switch (result) {
    case (Gtk::RESPONSE_OK): {
      cerr << "Opening file " << FileChooser.get_filename() << endl;
      Read (FileChooser.get_filename());
      FrameChanged();
      break;
    }
    case (Gtk::RESPONSE_CANCEL): {
      cerr << "Cancel.\n";
      break;
    }
  }
  
  FileChooser.hide();
}

void VisualClass::OnExport()
{
  cerr << "Export called.\n";
}

void VisualClass::ResetView()
{
  PathVis.View.Reset();
  double maxDim = max(max(Box[0], Box[1]), Box[2]);
  PathVis.View.SetDistance (1.2*maxDim);
  PathVis.Invalidate();
}

// bool VisualClass::on_delete_event()
// {
//   cerr << "delete event called.\n";
//   return true;
// }

void VisualClass::Quit()
{
  Gtk::Main::quit();
}

void VisualClass::FrameChanged()
{
  MakeFrame ((int)floor(FrameAdjust.get_value()));
  PathVis.Invalidate();
}


//////////
// Main.//
//////////

int main(int argc, char** argv)
{
  Gtk::Main kit(argc, argv);

  // Init gtkglextmm.
  Gtk::GL::init(argc, argv);

  if (argc < 2) {
    cerr << "Usage:\n  Visual myfile.h5\n";
    exit (1);
  }
  

  // Query OpenGL extension version.
  int major, minor;
  Gdk::GL::query_version(major, minor);
  std::cout << "OpenGL extension version - "
            << major << "." << minor << std::endl;

//   // Instantiate and run the application.
  VisualClass visual;
  visual.Read (argv[1]);
  kit.run(visual);

//   Vec3 r2 (0.99, 0.2, 0.3), r1(1.09, 0.3, 0.4), wall1, wall2;
//   BoxClass box;
//   box.Set(2.0, 2.0, 2.0);
//   if (box.BreakSegment(r1, r2, wall1, wall2)) {
//     cerr << "r1    = " << r1 << endl;
//     cerr << "r2    = " << r2 << endl;
//     cerr << "wall1 = " << wall1 << endl;
//     cerr << "wall2 = " << wall2 << endl;
//   }

  return 0;
}


void VisualClass::LineToggle()
{
  if (LinesButton.get_active())
    PathType = LINES;
  else
    PathType = TUBES;

  FrameChanged();
}


void VisualClass::MakePaths(int frame)
{
  int numPtcls  = PathArray.extent(1);
  int numSlices = PathArray.extent(2);
  Array<bool,1> used(numPtcls);
  used = false;

  for (int i=0; i<Paths.size(); i++)
    delete Paths[i];
  Paths.resize(0);

  vector<vector<int> > loopList;

  // Constuct list of permuting loops.  Ignore classical particles.
  for (int ptcl=0; ptcl<numPtcls; ptcl++) 
    if (!used(ptcl) && (PtclSpecies(ptcl).lambda!=0.0)) {
      vector<int> loop;
      int permPtcl = ptcl;
       do {
	 loop.push_back(permPtcl);
	 used(permPtcl) = true;
	 permPtcl = PermArray(frame, permPtcl);
       } while (permPtcl != ptcl);
      loopList.push_back(loop);
    }

  //  Paths.resize(loopList.size());
  cerr << loopList.size() << " loops.\n";
  for (int li=0; li<loopList.size(); li++) {
    vector<int> &loop = loopList[li];
    OnePath &path = (*new OnePath);
    path.Path.resize(numSlices*loop.size()+1);
    path.Color.resize(numSlices*loop.size()+1);
    int offset = 0;
    for (int pi=0; pi<loop.size(); pi++) {
      int ptcl = loop[pi];
      for (int slice=0; slice < numSlices; slice++) {
	path.Path[slice+offset][0] = PathArray(frame,ptcl,slice,0);
	path.Path[slice+offset][1] = PathArray(frame,ptcl,slice,1);
	path.Path[slice+offset][2] = PathArray(frame,ptcl,slice,2);
      }
      offset +=  numSlices;
    }
    // Close the path!!!
    path.Path[path.Path.size()-1][0] = PathArray(frame,loop[0],0,0);
    path.Path[path.Path.size()-1][1] = PathArray(frame,loop[0],0,1);
    path.Path[path.Path.size()-1][2] = PathArray(frame,loop[0],0,2);
    // Now, make sure the path doesn't have any discontinuities 
    // because having different period images.
    for (int slice=0; slice<path.Path.size()-1; slice++) 
      for (int dim=0; dim<3; dim++) {
	while ((path.Path[slice+1][dim] - path.Path[slice][dim]) >0.5*Box[dim])
	  path.Path[slice+1][dim] -= Box[dim];
	while ((path.Path[slice+1][dim] - path.Path[slice][dim])<-0.5*Box[dim])
	  path.Path[slice+1][dim] += Box[dim];
      }
    Paths.push_back (&path);
  }
  if (Wrap)
    for (int i=0; i<3; i++)
      Box.PutPathsInBox (Paths);

}



// Returns true if a segment walks across a box boundary.  In this
// case, it "breaks" the segment into two disjoint pieces, one from
// r1 to wall1 and one from wall2 to r2.  Note that to account for
// segments that cross more than one wall (i.e. a corner-crosser), we
// have to call this 3 times to ensure we get it right.
bool BoxClass::BreakSegment (Vec3 &r1, Vec3 &r2,
			     Vec3 &wall1, Vec3 &wall2)
{
  PutInBox (r1);
  PutInBox (r2);
  dVec disp = r2-r1;
  double eps = 1.0e-12;
  for (int dim=0; dim<3; dim++) 
    if (disp[dim]<-0.5*Box[dim]) { // path wraps in + direction
      double d1 = 0.5*Box[dim]-r1[dim];
      double d2 = r2[dim] + 0.5*Box[dim];

      double s1 = d2/(d1+d2);
      double s2 = 1.0-s1;
      Vec3 rshift = r2;
      rshift[dim] += Box[dim];
      wall1 = s1*r1 + s2*rshift;
      wall2 = wall1;
      wall2[dim] -= Box[dim];
      wall1[dim] -= eps;
      wall2[dim] += eps;
      return true;
    }
    else if (disp[dim]>0.5*Box[dim]) {
      double d2 = 0.5*Box[dim] - r2[dim];
      double d1 = r1[dim] + 0.5*Box[dim];
      double s1 = d2/(d1+d2);
      double s2 = 1.0-s1;
      Vec3 rshift = r2;
      rshift[dim] -= Box[dim];
      wall1 = s1*r1 + s2*rshift;
      wall2 = wall1;
      wall2[dim] += Box[dim];
      wall1[dim] += eps;
      wall2[dim] -= eps;
      return true;
    }
  return false;
}


bool BoxClass::BreakSegment (Vec3 &r1, Vec3 &r2,
			     Vec3 &wall1, Vec3 &wall2,
			     int dim)
{
  PutInBox (r1, dim);
  PutInBox (r2, dim);
  dVec disp = r2-r1;
  double eps = 1.0e-12;
  if (disp[dim]<-0.5*Box[dim]) { // path wraps in + direction
    double d1 = 0.5*Box[dim]-r1[dim];
    double d2 = r2[dim] + 0.5*Box[dim];
    
    double s1 = d2/(d1+d2);
    double s2 = 1.0-s1;
    Vec3 rshift = r2;
    rshift[dim] += Box[dim];
    wall1 = s1*r1 + s2*rshift;
    wall2 = wall1;
    wall2[dim] -= Box[dim];
    wall1[dim] -= eps;
    wall2[dim] += eps;
    return true;
  }
  else if (disp[dim]>0.5*Box[dim]) {
    double d2 = 0.5*Box[dim] - r2[dim];
    double d1 = r1[dim] + 0.5*Box[dim];
    double s1 = d2/(d1+d2);
    double s2 = 1.0-s1;
    Vec3 rshift = r2;
    rshift[dim] -= Box[dim];
    wall1 = s1*r1 + s2*rshift;
    wall2 = wall1;
    wall2[dim] += Box[dim];
    wall1[dim] += eps;
    wall2[dim] -= eps;
    return true;
  }
  return false;
}



void BoxClass::PutPathsInBox (vector<OnePath*>& inList)
{
  // To be perfectly correct, we need to repeat this process three
  // times to ensure that segments that cross the box corners get
  // properly broken up into two or three segments.
  for (int dim=0; dim<3; dim++) {
    vector<OnePath*> outList;
    outList.resize(0);
    for (int in=0; in<inList.size(); in++) {
      OnePath *newPath = new OnePath;
      OnePath *oldPath = inList[in];
      for (int slice=0; slice<(oldPath->Path.size()-1); slice++) {
	Vec3 r1, r2, wall1, wall2;
	r1 = oldPath->Path[slice];
	r2 = oldPath->Path[slice+1];
	if (BreakSegment (r1, r2, wall1, wall2, dim)) {
	  newPath->Path.push_back(r1);
	  newPath->Path.push_back(wall1);
	  outList.push_back(newPath);
	  newPath = new OnePath;
	  newPath->Path.push_back(wall2);
	} 
	else
	  newPath->Path.push_back(r1);
      }
      // Push last coordinate onto new path
      Vec3 rlast = oldPath->Path[oldPath->Path.size()-1];
      PutInBox(rlast, dim);
      newPath->Path.push_back(rlast);
      // Close the path
//       rlast = oldPath->Path[0];
//       PutInBox(rlast, dim);
//       newPath->Path.push_back(rlast);
      outList.push_back (newPath);
      delete oldPath;
    }
    inList.resize(outList.size());
    for (int i=0; i<outList.size(); i++)
      inList[i] = outList[i];
  }
}


void VisualClass::WrapToggle()
{
  Wrap = WrapButton.get_active();
  FrameChanged();
}
