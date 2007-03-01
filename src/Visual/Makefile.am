bin_PROGRAMS = pathvis++ mdvis++ wfvis++

LDADD = $(GTKMM_LIBS) $(GTKGLEXTMM_LIBS) $(COMMON_LIBS) $(GTS_LIBS) $(GTHREAD_LIBS)
AM_CXXFLAGS = $(GTKMM_CFLAGS) $(GTKGLEXTMM_CFLAGS) $(COMMON_CFLAGS)\
              $(GTS_CFLAGS) $(GTHREAD_CFLAGS) -DPKG_DATA_DIR=\"$(pkgdatadir)/\"
pathvis___SOURCES =       \
    BoxClass.cc         \
    BoxClass.h          \
    BoxObject.cc        \
    BoxObject.h         \
    ColorMap.cc         \
    ColorMap.h          \
    CoordObject.cc      \
    CoordObject.h       \
    DiskObject.cc       \
    DiskObject.h        \
    Export.cc           \
    Export.h            \
    GLObject.cc         \
    GLObject.h          \
    OnePath.h           \
    PathObject.cc       \
    PathObject.h        \
    PathVis.cc          \
    PathVis.h           \
    PlaneObject.cc      \
    PlaneObject.h       \
    SmoothClass.cc      \
    SmoothClass.h       \
    SphereObject.cc     \
    SphereObject.h      \
    ViewClass.cc        \
    Visual.cc           \
    Visual.h            \
    Isosurface.cc       \
    Isosurface.h        \
    trackball.h         \
    trackball.c    

mdvis___SOURCES =       \
    BoxClass.cc         \
    BoxClass.h          \
    BoxObject.cc        \
    BoxObject.h         \
    CoordObject.cc      \
    CoordObject.h       \
    DiskObject.cc       \
    DiskObject.h        \
    Isosurface.h        \
    Isosurface.cc       \
    MDExport.cc         \
    MDExport.h          \
    GLObject.cc         \
    GLObject.h          \
    MDVis.cc            \
    MDVis.h             \
    PathVis.cc          \
    PathVis.h           \
    SphereObject.cc     \
    SphereObject.h      \
    ViewClass.cc        \
    trackball.h         \
    trackball.c    

wfvis___SOURCES =       \
    BoxClass.cc         \
    BoxClass.h          \
    BoxObject.cc        \
    BoxObject.h         \
    ColorMap.h          \
    ColorMap.cc         \
    CoordObject.cc      \
    CoordObject.h       \
    DiskObject.cc       \
    DiskObject.h        \
    ElementData.cc      \
    ElementData.h       \
    Isosurface.h        \
    Isosurface.cc       \
    GLObject.cc         \
    GLObject.h          \
    WFVis.cc            \
    WFVis.h             \
    ParseCommand.cc     \
    ParseCommand.h      \
    PathVis.cc          \
    PathVis.h           \
    PlaneObject.cc      \
    PlaneObject.h       \
    SphereObject.cc     \
    SphereObject.h      \
    ViewClass.cc        \
    WFExport.cc         \
    WFExport.h          \
    trackball.h         \
    trackball.c    


dist_pkgdata_DATA = tubes.png lines.png straight.png smooth.png nowrap2.png \
                    wrap.png orthographic.png perspective.png player_play.png \
                    player_pause.png player_rev.png clipping.png              \
                    isoButton.png colormaps.in
