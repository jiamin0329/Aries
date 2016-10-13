
message(STATUS "cmake dualgrid")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")
include_directories(../common)
include_directories(../procdata)

set(DGRD_SRC DualGrid.cpp DGEdge.cpp DGVertex.cpp DGPoint.cpp)

add_library(libdgrd ${DGRD_SRC})
set_target_properties(libdgrd PROPERTIES OUTPUT_NAME "dgrd")
