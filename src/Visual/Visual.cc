#include "Visual.h"

void VisualClass::MakeFrame(int frame)
{
  int numSpecies = Species.size();

  for (int i=0; i<PathVis.Objects.size(); i++)
    delete PathVis.Objects[i];
  
  PathVis.Objects.resize(0);

  int numPtcls = Paths.extent(1);
  int numSlices = Paths.extent(2);


  Array<Vec3,1> onePath(numSlices);
  for (int si=0; si<numSpecies; si++) {
    for (int ptcl=Species(si).FirstParticle; 
	 ptcl<=Species(si).LastParticle; ptcl++) {
      if (Species(si).lambda != 0.0) {
	for (int slice=0; slice<numSlices; slice++) {
	  onePath(slice)[0] = Paths(frame,ptcl,slice,0);
	  onePath(slice)[1] = Paths(frame,ptcl,slice,1);
	  onePath(slice)[2] = Paths(frame,ptcl,slice,2);
	}
	PathObject* pathObj = new PathObject;
	if (si == 0)
	  pathObj->SetColor (0.3, 0.3, 1.0);
	else
	  pathObj->SetColor (0.0, 1.0, 0.0);
	if (PathType == TUBES)
	  pathObj->TubesSet (onePath);
	else
	  pathObj->LinesSet (onePath);
	PathVis.Objects.push_back(pathObj);
      }
      else {
	Vec3 pos;
	pos[0] = Paths(frame, ptcl, 0, 0);
	pos[1] = Paths(frame, ptcl, 0, 1);
	pos[2] = Paths(frame, ptcl, 0, 2);

	SphereObject* sphere = new SphereObject;
	sphere->SetPos (pos);
	if (ptcl != 30)
	  sphere->SetColor (Vec3(1.0, 0.0, 0.0));
	else
	  sphere->SetColor (Vec3(1.0, 0.0, 1.0));
	PathVis.Objects.push_back(sphere);
      }
    }
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
  assert(in.ReadVar ("Path", Paths));
  PutInBox();

  FrameAdjust.set_upper(Paths.extent(0)-1);
  
  in.CloseSection();
  in.CloseFile();
  MakeFrame (0);
}

void VisualClass::PutInBox()
{
  for (int frame=0; frame<Paths.extent(0); frame++) 
    for (int ptcl=0; ptcl<Paths.extent(1); ptcl++) 
      for (int dim=0; dim<3; dim++) {
	while (Paths(frame,ptcl,0,dim) > 0.5*Box[dim])
	  for (int slice=0; slice<Paths.extent(2); slice++)
	    Paths(frame,ptcl,slice,dim) -= Box[dim];
	while (Paths(frame,ptcl,0,dim) < -0.5*Box[dim])
	  for (int slice=0; slice<Paths.extent(2); slice++)
	    Paths(frame,ptcl,slice,dim) += Box[dim];
      }
	  
}


#include "tubes.xpm"
#include "lines.xpm"

VisualClass::VisualClass()
  : m_VBox(false, 0), m_ButtonQuit("Quit"), 
    FrameAdjust (0.0, 0.0, 0.0),
    PathType (LINES),
    TubesImage("tubes.png"), LinesImage("lines.png")
{
  // Top-level window.
  set_title("VisualClass");

  // Get automatically redrawn if any of their children changed allocation.
  set_reallocate_redraws(true);
  add(m_VBox);

  // VisualClass OpenGL scene.
  PathVis.set_size_request(600, 600);


  // VisualClass quit button.

  m_ButtonQuit.signal_clicked().connect(
    sigc::mem_fun(*this, &VisualClass::on_button_quit_clicked));
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
  Tools.append (LinesButton);
  Tools.append (TubesButton);


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
		sigc::mem_fun(*this, &VisualClass::on_delete_event));
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
  cerr << "Open called\n";
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

void VisualClass::on_delete_event()
{
  cerr << "delete event called.\n";
}

void VisualClass::on_button_quit_clicked()
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
