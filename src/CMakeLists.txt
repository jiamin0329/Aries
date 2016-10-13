
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin")
set(LIBRARY_OUTPUT_PATH    "${PROJECT_SOURCE_DIR}/lib")

#set link directory
link_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories(./common)
include_directories(./procdata)

add_executable(ARIES ${PROJECT_SOURCE_DIR}/src/ARIES_main.cpp)
target_link_libraries(ARIES comm dgrd mat geom pdat grid)
