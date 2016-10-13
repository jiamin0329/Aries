




message(STATUS "cmake output")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

set(OPUT_SRC Output.cpp)

add_library(liboput ${OPUT_SRC})
set_target_properties(liboput PROPERTIES OUTPUT_NAME "oput")
