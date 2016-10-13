




message(STATUS "cmake material")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")
include_directories(../common)
include_directories(../procdata)

set(MAT_SRC ViscosityModel.cpp ConductivityModel.cpp FluidModel.cpp)

add_library(libmat ${MAT_SRC})
set_target_properties(libmat PROPERTIES OUTPUT_NAME "mat")
