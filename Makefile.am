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
