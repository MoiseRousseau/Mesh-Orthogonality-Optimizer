cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(OrthOpt LANGUAGES CXX VERSION 1.0 DESCRIPTION "Mesh Orthogonality optimizer")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math")

#------------------
# Define target 
#------------------
include_directories(src)
add_executable (OrthOpt 
    src/main.cpp
    src/Connection.cpp
    src/Mesh.cpp
    src/OrthOpt.cpp
    src/IO.cpp
)
#set_target_properties(OrthOpt PROPERTIES VERSION ${PROJECT_VERSION})

# Search Eigen3
message("Searching for Eigen3 library")
find_package(
    Eigen3 REQUIRED
    HINTS ${EIGEN3_PATH}
)
message("Eigen3 library found in: ${EIGEN3_INCLUDE_DIR}")

#include external dependencies (LBFGS++)
set( LBGFSpp "${CMAKE_SOURCE_DIR}/external" )
message("${LBGFSpp}")

target_include_directories ( OrthOpt PRIVATE ${EIGEN3_INCLUDE_DIR} ${LBGFSpp} )

#------------------
# Parallelization using OpenMP
#------------------
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(OrthOpt OpenMP::OpenMP_CXX)
endif()


# install(TARGETS qcm DESTINATION ${PYTHON_INCLUDE_DIRS}/site-packages)
# install(TARGETS qcm DESTINATION ${CMAKE_SOURCE_DIR}/pyqcm)
install(TARGETS OrthOpt DESTINATION OrthOpt)


