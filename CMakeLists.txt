cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(OrthOpt LANGUAGES CXX VERSION 1.0 DESCRIPTION "Mesh Orthogonality optimizer")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -ftree-vectorize -ffast-math")

#------------------
# Define target 
#------------------
include_directories(src)
add_executable (helloDemo 
    demo.cxx 
    demo_b.cxx
)
set_target_properties(qcm PROPERTIES VERSION ${PROJECT_VERSION} PREFIX "" SUFFIX ".so")

set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()
include_directories ( "${EIGEN3_INCLUDE_DIR}" )

#------------------
# Link dependencies
#------------------
target_link_libraries(qcm ${PYTHON_LIBRARIES} ${CUBA_LIBRARY} ${LAPACK_LIBRARIES})
target_include_directories(qcm PUBLIC src_qcm src_ed src_python src_util ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR}/numpy ${CUBA_DIR})


#------------------
# Parallelization using OpenMP
#------------------
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(OrthOpt OpenMP::OpenMP_CXX)
endif()


# install(TARGETS qcm DESTINATION ${PYTHON_INCLUDE_DIRS}/site-packages)
# install(TARGETS qcm DESTINATION ${CMAKE_SOURCE_DIR}/pyqcm)
install(TARGETS qcm DESTINATION pyqcm)

add_test(OrthOpt)

