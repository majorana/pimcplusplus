#include "Visual.h"

void VisualClass::Read(string fileName)
{
  IOSectionClass in;
  assert(in.OpenFile (fileName));

  Array<double,1> box;
  assert (in.OpenSection("System"));
  assert (in.ReadVar ("Box", box));
  cerr << "Box = " << box << endl;
  BoxObject *boxObject = new BoxObject;
  boxObject->Set (box(0), box(1), box(2));
  boxObject->SetColor (0.5, 0.5, 1.0);
  double maxDim = max(max(box(0), box(1)), box(2));
  PathVis.Objects.push_back(boxObject);
  //  PathVis.View.SetDistance (1.2*maxDim);
  PathVis.View.SetDistance (0.2*maxDim);

  int numSpecies = in.CountSections ("Species");  
  in.CloseSection (); // "System"

  assert(in.OpenSection("Observables"));
  assert(in.OpenSection("PathDump"));
  Array<double,4> paths;
  assert(in.ReadVar ("Path", paths));
  
  int numPtcls = paths.extent(1);
  int numSlices = paths.extent(2);

  cerr << "numPtcls = " << numPtcls << endl;
  cerr << "numSlices = " << numSlices << endl;

  Array<Vec3,1> onePath(numSlices);
  for (int ptcl=0; ptcl<numPtcls; ptcl++) {
    for (int slice=0; slice<numSlices; slice++) {
//       fprintf (stderr, "[ %8.5f %8.5f %8.5f ]\n",
// 	       paths(0,ptcl,slice,0),
// 	       paths(0,ptcl,slice,1),
// 	       paths(0,ptcl,slice,2));
      onePath(slice)[0] = paths(500,ptcl,slice,0);
      onePath(slice)[1] = paths(500,ptcl,slice,1);
      onePath(slice)[2] = paths(500,ptcl,slice,2);
    }
    PathObject* pathObj = new PathObject;
    pathObj->Set (onePath);
    pathObj->SetColor (0.3, 0.3, 1.0);
    PathVis.Objects.push_back(pathObj);
  }

  in.CloseSection();

  in.CloseFile();
}

VisualClass::VisualClass()
  : m_VBox(false, 0), m_ButtonQuit("Quit")
{
  // Top-level window.
  set_title("VisualClass");

  // Get automatically redrawn if any of their children changed allocation.
  set_reallocate_redraws(true);
  add(m_VBox);

  // VisualClass OpenGL scene.
  PathVis.set_size_request(400, 400);

  m_VBox.pack_start(PathVis);
  // VisualClass quit button.

  m_ButtonQuit.signal_clicked().connect(
    sigc::mem_fun(*this, &VisualClass::on_button_quit_clicked));
  m_VBox.pack_start(m_ButtonQuit, Gtk::PACK_SHRINK, 0);

  // Show window.
  show_all();
}

VisualClass::~VisualClass()
{}

void VisualClass::on_button_quit_clicked()
{
  Gtk::Main::quit();
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

  Array<Vec3, 1> path(5);
  path(0) = Vec3(-0.5,  0.5, -0.5);
  path(1) = Vec3( 0.5,  0.5, -0.5);
  path(2) = Vec3( 0.5, -0.5, -0.5);
  path(3) = Vec3(-0.5, -0.5, -0.5);
  path(4) = Vec3(-0.5,  0.5, -0.5);

  PathObject *p1 = new PathObject();
  p1->SetColor (0.0, 0.0, 1.0);
  p1->Set (path);

  visual.PathVis.Objects.push_back(p1);

  // visual.PathVis.AddPath (path);
  for (int i=0; i<5; i++)
    path(i) += Vec3(0.0, 0.0, 1.0);
  //  visual.PathVis.AddPath (path);
  
  PathObject *p2 = new PathObject();
  p2->SetColor (1.0, 0.0, 0.0);
  p2->Set (path);
  visual.PathVis.Objects.push_back(p2);
  
  BoxObject *box = new BoxObject;
  box->Set (2.0, 1.0, 0.5);
  visual.PathVis.Objects.push_back(box);


  kit.run(visual);

  return 0;
}
