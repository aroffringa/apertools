cmake_minimum_required(VERSION 2.8.8)

project(apertools)

add_definitions(-DAOPROJECT)

if(PORTABLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wvla -Wall -DNDEBUG -std=c++11")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wvla -Wall -DNDEBUG -march=native -std=c++11")
endif(PORTABLE)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits scimath scimath_f)

find_package(CFITSIO REQUIRED)

find_library(PTHREAD_LIB pthread REQUIRED)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${CFITSIO_INCLUDE_DIR})

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

add_executable(apbeam apbeam.cpp fitsreader.cpp fitswriter.cpp fitsiochecker.cpp image.cpp)
target_link_libraries(apbeam ${CASACORE_LIBRARIES} ${CFITSIO_LIBRARY})

add_executable(applybeam applybeam.cpp fitsreader.cpp fitswriter.cpp fitsiochecker.cpp)
target_link_libraries(applybeam ${CASACORE_LIBRARIES} ${CFITSIO_LIBRARY})

message(STATUS "Flags passed to C++ compiler: " ${CMAKE_CXX_FLAGS})
