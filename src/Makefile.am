SUBDIRS = Common
bin_PROGRAMS = pimc++2D pimc++3D TestPerm TestEwald

pimc++2D_CFLAGS = -DTWO_D
pimc++3D_CFLAGS = -DTHREE_D
pimc++2D_SOURCES =                    \
  Main.cc                             \
  ObservableEnergy.cc                 \
  BisectionClass.cc                   \
  PIMCClass.cc                        \
  MetaMoves.cc                        \
  BlockMove.cc                        \
  ObservableClass.cc                  \
  Common/Splines/CubicSpline.cc       \
  Common/Splines/MyTricubicSpline.cc  \
  Common/Splines/Grid.cc              \
  Common/Splines/QuinticSpline.cc     \
  Common/Splines/QuinticSplines.cc    \
  SpeciesClass.cc                     \
  Common.cc                           \
  PermuteTableClass.cc                \
  BisectionMoveClass.cc               \
  MoveClass.cc                        \
  ActionClass.cc                      \
  PathDataClass.cc                    \
  CommunicatorClass.cc                \
  PathClass.cc                        \
  MirroredArrayClass.cc               \
  WrapClass.cc                        \
  Common/MPI/Communication.cc         \
  Common/IO/InputOutput.cc            \
  Common/IO/InputOutputHDF5.cc        \
  Common/IO/InputFile.cc              \
  Common/IO/InputOutputASCII.cc       \
  Common/IO/InputOutputXML.cc         \
  Common/PairAction/PAcoulombBCFit.cc \
  Common/PairAction/PAclassicalFit.cc \
  Common/PairAction/PAzeroFit.cc      \
  Common/Splines/BicubicSpline.cc     \
  Common/PH/Potential.cc              \
  Common/PH/QuinticPH.cc              \
  Common/PH/CoulombPot.cc             \
  Common/PH/ScreenedPot.cc            \
  Common/PH/SplinePot.cc              \
  Common/PH/HeAzizPot.cc              \
  Common/PairAction/DavidPAClass.cc   \
  Common/Ewald/OptimizedBreakup.cc    \
  Common/MatrixOps/MatrixOps.cc       \
  Common/Integration/GKIntegration.cc \
  MirroredClass.cc

pimc++3D_SOURCES =                    \
  Main.cc                             \
  ObservableEnergy.cc                 \
  BisectionClass.cc                   \
  PIMCClass.cc                        \
  MetaMoves.cc                        \
  BlockMove.cc                        \
  ObservableClass.cc                  \
  Common/Splines/CubicSpline.cc       \
  Common/Splines/MyTricubicSpline.cc  \
  Common/Splines/Grid.cc              \
  Common/Splines/QuinticSpline.cc     \
  Common/Splines/QuinticSplines.cc    \
  SpeciesClass.cc                     \
  Common.cc                           \
  PermuteTableClass.cc                \
  BisectionMoveClass.cc               \
  MoveClass.cc                        \
  ActionClass.cc                      \
  PathDataClass.cc                    \
  CommunicatorClass.cc                \
  PathClass.cc                        \
  MirroredArrayClass.cc               \
  WrapClass.cc                        \
  Common/MPI/Communication.cc         \
  Common/IO/InputOutput.cc            \
  Common/IO/InputOutputHDF5.cc        \
  Common/IO/InputFile.cc              \
  Common/IO/InputOutputASCII.cc       \
  Common/IO/InputOutputXML.cc         \
  Common/PairAction/PAcoulombBCFit.cc \
  Common/PairAction/PAclassicalFit.cc \
  Common/PairAction/PAzeroFit.cc      \
  Common/Splines/BicubicSpline.cc     \
  Common/PH/Potential.cc              \
  Common/PH/QuinticPH.cc              \
  Common/PH/CoulombPot.cc             \
  Common/PH/ScreenedPot.cc            \
  Common/PH/SplinePot.cc              \
  Common/PH/HeAzizPot.cc              \
  Common/PairAction/DavidPAClass.cc   \
  Common/Ewald/OptimizedBreakup.cc    \
  Common/MatrixOps/MatrixOps.cc       \
  Common/Integration/GKIntegration.cc \
  MirroredClass.cc


TestPerm_SOURCES =                    \
  TestPermutation.cc                  \
  BisectionClass.cc                   \
  BlockMove.cc                        \
  MetaMoves.cc                        \
  PIMCClass.cc                        \
  ObservableClass.cc                  \
  ObservableEnergy.cc                 \
  Common/Splines/CubicSpline.cc       \
  Common/Splines/MyTricubicSpline.cc  \
  Common/Splines/Grid.cc              \
  Common/Splines/QuinticSpline.cc     \
  Common/Splines/QuinticSplines.cc    \
  SpeciesClass.cc                     \
  Common.cc                           \
  PermuteTableClass.cc                \
  BisectionMoveClass.cc               \
  MoveClass.cc                        \
  ActionClass.cc                      \
  PathDataClass.cc                    \
  CommunicatorClass.cc                \
  PathClass.cc                        \
  MirroredArrayClass.cc               \
  WrapClass.cc                        \
  Common/MPI/Communication.cc         \
  Common/IO/InputOutput.cc            \
  Common/IO/InputOutputHDF5.cc        \
  Common/IO/InputFile.cc              \
  Common/IO/InputOutputASCII.cc       \
  Common/IO/InputOutputXML.cc         \
  Common/PairAction/PAcoulombBCFit.cc \
  Common/PairAction/PAclassicalFit.cc \
  Common/PairAction/PAzeroFit.cc      \
  Common/Splines/BicubicSpline.cc     \
  Common/PH/Potential.cc              \
  Common/PH/QuinticPH.cc              \
  Common/PH/CoulombPot.cc             \
  Common/PH/ScreenedPot.cc            \
  Common/PH/SplinePot.cc              \
  Common/PH/HeAzizPot.cc              \
  Common/PairAction/DavidPAClass.cc   \
  Common/Ewald/OptimizedBreakup.cc    \
  Common/MatrixOps/MatrixOps.cc       \
  Common/Integration/GKIntegration.cc \
  MirroredClass.cc

TestEwald_SOURCES =                   \
  TestEwald.cc                        \
  PathClass.cc                        \
  MirroredClass.cc                    \
  SpeciesClass.cc                     \
  ActionClass.cc                      \
  PathDataClass.cc                    \
  Common/MPI/Communication.cc         \
  Common/IO/InputOutput.cc            \
  Common/IO/InputOutputHDF5.cc        \
  Common/IO/InputFile.cc              \
  Common/IO/InputOutputASCII.cc       \
  Common/IO/InputOutputXML.cc         \
  Common/PairAction/PAcoulombBCFit.cc \
  Common/PairAction/PAclassicalFit.cc \
  Common/PairAction/PAzeroFit.cc      \
  Common/PairAction/DavidPAClass.cc   \
  Common/Splines/BicubicSpline.cc     \
  Common/PH/Potential.cc              \
  Common/PH/QuinticPH.cc              \
  Common/PH/CoulombPot.cc             \
  Common/PH/ScreenedPot.cc            \
  Common/PH/SplinePot.cc              \
  Common/PH/HeAzizPot.cc              \
  Common/Splines/CubicSpline.cc       \
  Common/Splines/MyTricubicSpline.cc  \
  Common/Splines/Grid.cc              \
  Common/Splines/QuinticSpline.cc     \
  Common/Ewald/OptimizedBreakup.cc    \
  Common/MatrixOps/MatrixOps.cc       \
  Common/Integration/GKIntegration.cc \
  Common/Splines/QuinticSplines.cc
