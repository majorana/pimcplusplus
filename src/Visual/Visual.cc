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
  Box[0] = box(0); Box[1] = box(1); Box[2] = box(2);
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
    FileChooser ("Choose an output file")
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
  LinesButton.set_label("Lines");
  TubesButton.set_label("Tubes");
  LinesButton.signal_toggled().connect
    (sigc::mem_fun(*this, &VisualClass::LineToggle));
  TubesButton.set_group (group);
  TubesButton.set_icon_widget (TubesImage);
  LinesButton.set_icon_widget (LinesImage);
  StraightButton.set_label("Straight");
  StraightButton.set_icon_widget(StraightImage);
  SmoothButton.set_label("Smooth");
  SmoothButton.set_icon_widget(SmoothImage);
  group = StraightButton.get_group();
  SmoothButton.set_group(group);
  Tools.append (LinesButton);
  Tools.append (TubesButton);
  Tools.append (ToolSep);
  Tools.append (StraightButton);
  Tools.append (SmoothButton);

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

  // Instantiate and run the application.
  VisualClass visual;
  visual.Read (argv[1]);
  kit.run(visual);

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
    path.Path.resize(numSlices*loop.size());
    path.Color.resize(numSlices*loop.size());
    int offset = 0;
    for (int pi=0; pi<loop.size(); pi++) {
      int ptcl = loop[pi];
      for (int slice=0; slice < numSlices; slice++) {
	path.Path(slice+offset)[0] = PathArray(frame,ptcl,slice,0);
	path.Path(slice+offset)[1] = PathArray(frame,ptcl,slice,1);
	path.Path(slice+offset)[2] = PathArray(frame,ptcl,slice,2);
      }
      offset +=  numSlices;
    }
    // Now, make sure the path doesn't have any discontinuities 
    // because having different period images.
    for (int slice=0; slice<path.Path.size()-1; slice++) 
      for (int dim=0; dim<3; dim++) {
	while ((path.Path(slice+1)[dim] - path.Path(slice)[dim]) >0.5*Box[dim])
	  path.Path(slice+1)[dim] -= Box[dim];
	while ((path.Path(slice+1)[dim] - path.Path(slice)[dim])<-0.5*Box[dim])
	  path.Path(slice+1)[dim] += Box[dim];
      }
    Paths.push_back (&path);
  }
  used = false;
  

}
