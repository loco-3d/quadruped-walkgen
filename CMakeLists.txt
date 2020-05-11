CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Set up project properties
SET(PROJECT_NAMESPACE loco-3d)
SET(PROJECT_NAME quadruped-walkgen)
SET(PROJECT_DESCRIPTION "MPC Controller using Crocoddyl")
SET(PROJECT_URL https://gitlab.laas.fr/${PROJECT_NAMESPACE}/${PROJECT_NAME})

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

# Disable -Werror on Unix.
SET(CXX_DISABLE_WERROR True)

ADD_PROJECT_DEPENDENCY(crocoddyl REQUIRED)

OPTION(BUILD_WITH_MULTITHREADS "Build the library with the Multithreading support (required OpenMP)" ON)
IF(BUILD_WITH_MULTITHREADS)
  SET(BUILD_WITH_NTHREADS "4" CACHE STRING "Number of threads")
  string(REGEX MATCH "^[0-9]+$" BUILD_WITH_NTHREADS  ${BUILD_WITH_NTHREADS})
  IF(NOT BUILD_WITH_NTHREADS MATCHES "^[0-9]+$")
    SET(BUILD_WITH_NTHREADS 4)
    MESSAGE(WARNING "the number of threads have to be an interger value, set to ${BUILD_WITH_NTHREADS}")
  ENDIF()
ENDIF()

FINDPYTHON(3.6 REQUIRED)
STRING(REGEX REPLACE "-" "_" PYTHON_DIR ${PROJECT_NAME})

SET(BOOST_REQUIERED_COMPONENTS python filesystem system)
SET(BOOST_COMPONENTS ${BOOST_REQUIERED_COMPONENTS})
SEARCH_FOR_BOOST()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include/${PROJECT_NAME})
ADD_SUBDIRECTORY(python)
