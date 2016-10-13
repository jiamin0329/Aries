#------------------------------------------------------------------------------
# - Set up FVM grid related source files
#
message(STATUS "cmake fvm grid")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

set(GRID_SRC
  Grid.cpp
  GridVertexMPI.cpp
  GridLine.cpp
  GridTriangle.cpp
  GridQuad.cpp
  GridHexa.cpp
  GridPrism.cpp
  GridTetra.cpp
  GridPyramid.cpp)

add_library(libgrid ${GRID_SRC})
set_target_properties(libgrid PROPERTIES OUTPUT_NAME "grid")
