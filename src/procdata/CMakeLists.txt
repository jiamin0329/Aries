
message(STATUS "cmake processor data")

set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")
include_directories(../common)

set(PROC_SRC ProcData.cpp)

add_library(libpdat ${PROC_SRC})
set_target_properties(libpdat PROPERTIES OUTPUT_NAME "pdat")
