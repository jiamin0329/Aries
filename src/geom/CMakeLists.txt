
message(STATUS "cmake geom")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

include_directories(../common)
include_directories(../dualgrid)
include_directories(../procdata)

set(GEOM_SRC MeshData.cpp )

add_library(libgeom ${GEOM_SRC})
set_target_properties(libgeom PROPERTIES OUTPUT_NAME "geom")
