########################################################
#                                                      #
# CMake File for the Project                           #
#                                                      #
# Top level CMakeLists.txt file                        #
#                                                      #
# Author: Jiamin Xu                                    #
# Date:   27th Aug, 2014                               #
# Email:  jiamin.xu@live.com                           #
#                                                      #
########################################################

#------------------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(STATUS "ARIES requires an out-of-source build.")
  message(FATAL_ERROR "in-source build detected")
endif()

#------------------------------------------------------------------------------
# - Define CMake requirements and override make rules as needed
#
# minimum version required
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# project name
project(ARIES)

# version number
set(ARIES_VERSION 0.1)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
message(STATUS ${CMAKE_BUILD_TYPE})
message(STATUS ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS ${CMAKE_CXX_FLAGS_DEBUG})

add_subdirectory("${PROJECT_SOURCE_DIR}/src/common")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/grid")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/procdata")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/geom")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/dualgrid")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/material")

#add_subdirectory("${PROJECT_SOURCE_DIR}/src/output")
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
