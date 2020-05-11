CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Set up project properties
SET(PROJECT_NAMESPACE loco-3d)
SET(PROJECT_NAME quadruped-walkgen)
SET(PROJECT_DESCRIPTION "MPC Controller using Crocoddyl")
SET(PROJECT_URL https://gitlab.laas.fr/loco-3d/quadruped-walkgen)

SET(CMAKE_VERBOSE_MAKEFILE True)

# Check if the submodule cmake have been initialized
IF(NOT EXISTS "${CMAKE_SOURCE_DIR}/cmake/base.cmake")
  MESSAGE(FATAL_ERROR "\nPlease run the following command first:\ngit submodule update --init\n")
ENDIF()

# Include important cmake modules
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/python.cmake)

COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
#set (CMAKE_CXX_STANDARD 11)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})
#ADD_DEFINITIONS(-DWITH_CODEGEN)

# Set a default build type to 'Release' if none was specified
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
ENDIF()

# Disable -Werror on Unix.
SET(CXX_DISABLE_WERROR True)

ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.0.5")
ADD_REQUIRED_DEPENDENCY("pinocchio >= 2.2.1")
ADD_REQUIRED_DEPENDENCY("eigenpy")
ADD_REQUIRED_DEPENDENCY("example-robot-data >= 3.1.0")
ADD_REQUIRED_DEPENDENCY("crocoddyl >= 0.9")

link_directories(/local/installation/lib)

OPTION(BUILD_WITH_MULTITHREADS "Build the library with the Multithreading support (required OpenMP)" ON)
IF(BUILD_WITH_MULTITHREADS)
  SET(BUILD_WITH_NTHREADS "4" CACHE STRING "Number of threads")
  string(REGEX MATCH "^[0-9]+$" BUILD_WITH_NTHREADS  ${BUILD_WITH_NTHREADS})
  IF(NOT BUILD_WITH_NTHREADS MATCHES "^[0-9]+$")
    SET(BUILD_WITH_NTHREADS 4)
    MESSAGE(WARNING "the number of threads have to be an interger value, set to ${BUILD_WITH_NTHREADS}")
  ENDIF()
ENDIF()

FINDPYTHON(3.6 EXACT REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${CROCODDYL_INCLUDE_DIRS} ${PINOCCHIO_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(bindings)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

SET(BOOST_REQUIERED_COMPONENTS python filesystem system)
SET(BOOST_COMPONENTS ${BOOST_REQUIERED_COMPONENTS})
SEARCH_FOR_BOOST()


ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include/${PROJECT_NAME})
ADD_SUBDIRECTORY(bindings)
