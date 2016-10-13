##**********************************************************
## CMakeList file for Aries Solver
##**********************************************************
message(STATUS "cmake common")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

set(COMM_SRC
  AriesManager.cpp 
  AriesMPI.cpp
  Logger.cpp
  StartupShutdownManager.cpp
  Utilities.cpp
  Clock.cpp)

add_library(libcomm ${COMM_SRC})
set_target_properties(libcomm PROPERTIES OUTPUT_NAME "comm")
