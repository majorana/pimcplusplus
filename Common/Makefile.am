SUBDIRS =                  \
          MPI              \
	  IO               \
          Splines          \
          SpecialFunctions \
          Random           \
          Integration      \
          MatrixOps        \
          Fitting          \
          DFT              \
          PH               \
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

if PARALLEL
  lib_LTLIBRARIES = libmpicommon.la 
  pkgconfig_DATA  = mpicommon.pc
else
  lib_LTLIBRARIES = libcommon.la
  pkgconfig_DATA  = common.pc
endif 
  EXTRA_DIST      = common.pc.in mpicommon.pc.in   \
                    turing.conf turing_serial.conf \
                    ice.conf                       \
                    pv.conf pv_serial.conf         \
                    maranatha.conf                 \
                    tungsten.conf                  \
                    cobalt.conf                  

libcommon_la_SOURCES =
libcommon_la_LIBADD = IO/libIO.la                    \
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

libmpicommon_la_SOURCES =
libmpicommon_la_LIBADD = IO/libIO.la                    \
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
include_HEADERS = Blitz.h config.h nan.h

pkgconfigdir = $(libdir)/pkgconfig/

install-data-local:
	$(INSTALL_DATA) $(pkgconfig_DATA) $(pkgconfigdir)


