SUBDIRS =                  \
	  IO               \
          Splines          \
          SpecialFunctions \
          Random           \
          Integration      \
          MatrixOps        \
          Fitting          \
          DFT              \
          PH               \
          MPI              \
          Distributed      \
          Ewald            \
          FFT              \
          PairAction       \
          PlaneWavePHDFT   \
          Atom           

AUX_DIST              = $(ac_aux_dir)/config.guess \
                        $(ac_aux_dir)/config.sub \
                        $(ac_aux_dir)/install-sh \
                        $(ac_aux_dir)/ltconfig \
                        $(ac_aux_dir)/ltmain.sh \
                        $(ac_aux_dir)/mdate-sh \
                        $(ac_aux_dir)/missing \
                        $(ac_aux_dir)/mkinstalldirs

lib_LTLIBRARIES = libcommon.la
libcommon_la_SOURCES =
libcommon_la_LIBADD = IO/libciao.la                  \
                      Splines/libspline.la           \
                      SpecialFunctions/libspecfun.la \
                      Integration/libquad.la         \
                      MatrixOps/libmatops.la         \
                      Fitting/libfitting.la          \
                      DFT/libdft.la                  \
                      PH/libPH.la                    \
                      MPI/libcomm.la                 \
                      Distributed/libdistrib.la      \
                      Ewald/libbreakup.la            \
                      FFT/libfft.la                  \
                      PairAction/libpairact.la       \
                      PlaneWavePHDFT/libpw.la        \
                      Atom/libatom.la 

includedir=$(prefix)/include/Common
include_HEADERS = Blitz.h config.h

pkgconfigdir = $(libdir)/pkgconfig/
pkgconfig_DATA = common.pc

install-data-local:
	$(INSTALL_DATA) common.pc $(pkgconfigdir)

EXTRA_DIST = common.pc.in turing.conf
