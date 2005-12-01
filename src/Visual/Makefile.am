bin_PROGRAMS = pathvis++ mdvis++

LDADD = $(GTKMM_LIBS) $(GTKGLEXTMM_LIBS) $(COMMON_LIBS) $(GTS_LIBS)
AM_CXXFLAGS = $(GTKMM_CFLAGS) $(GTKGLEXTMM_CFLAGS) $(COMMON_CFLAGS)\
              $(GTS_CFLAGS) -DPKG_DATA_DIR=\"$(pkgdatadir)/\"
pathvis___SOURCES =       \
    BoxClass.cc         \
    BoxClass.h          \
    BoxObject.cc        \
    BoxObject.h         \
    Export.cc           \
    Export.h            \
    GLObject.cc         \
    GLObject.h          \
    OnePath.h           \
    PathObject.cc       \
    PathObject.h        \
    PathVis.cc          \
    PathVis.h           \
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


dist_pkgdata_DATA = tubes.png lines.png straight.png smooth.png nowrap2.png \
                    wrap.png orthographic.png perspective.png player_play.png \
                    player_pause.png player_rev.png
