CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Set up project properties
SET(PROJECT_NAMESPACE loco-3d)
SET(PROJECT_NAME quadruped-walkgen)
SET(PROJECT_DESCRIPTION "MPC Controller using Crocoddyl")
SET(PROJECT_URL "https://github.com/${PROJECT_NAMESPACE}/${PROJECT_NAME}")

OPTION(SUFFIX_SO_VERSION "Suffix library name with its version" ON)

SET(CMAKE_VERBOSE_MAKEFILE TRUE)
SET(CUSTOM_HEADER_DIR ${PROJECT_NAME})
SET(CXX_DISABLE_WERROR TRUE)
SET(PROJECT_USE_CMAKE_EXPORT TRUE)

# Check if the submodule cmake have been initialized
IF(NOT EXISTS "${CMAKE_SOURCE_DIR}/cmake/base.cmake")
  MESSAGE(FATAL_ERROR "\nPlease run the following command first:\ngit submodule update --init\n")
ENDIF()

# Include important cmake modules
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/python.cmake)

COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})
#ADD_DEFINITIONS(-DWITH_CODEGEN)

ADD_PROJECT_DEPENDENCY(crocoddyl REQUIRED)
ADD_PROJECT_DEPENDENCY(example-robot-data)
FIND_PACKAGE(Boost REQUIRED COMPONENTS filesystem system)
FINDPYTHON(REQUIRED)
SEARCH_FOR_BOOST_PYTHON(REQUIRED)
STRING(REGEX REPLACE "-" "_" PYTHON_DIR ${PROJECT_NAME})

OPTION(BUILD_WITH_MULTITHREADS "Build the library with the Multithreading support (required OpenMP)" ON)
IF(BUILD_WITH_MULTITHREADS)
  SET(BUILD_WITH_NTHREADS "4" CACHE STRING "Number of threads")
  string(REGEX MATCH "^[0-9]+$" BUILD_WITH_NTHREADS  ${BUILD_WITH_NTHREADS})
  IF(NOT BUILD_WITH_NTHREADS MATCHES "^[0-9]+$")
    SET(BUILD_WITH_NTHREADS 4)
    MESSAGE(WARNING "the number of threads have to be an interger value, set to ${BUILD_WITH_NTHREADS}")
  ENDIF()
ENDIF()

SET(${PROJECT_NAME}_HEADERS
  include/${CUSTOM_HEADER_DIR}/quadruped_augmented.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_augmented.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped_augmented_time.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_augmented_time.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped_nl.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_nl.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped_step.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_step.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped_step_period.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_step_period.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped_step_time.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_step_time.hxx
  include/${CUSTOM_HEADER_DIR}/quadruped_time.hpp
  include/${CUSTOM_HEADER_DIR}/quadruped_time.hxx
  )

SET(${PROJECT_NAME}_SOURCES
  src/quadruped.cpp
  src/quadruped_nl.cpp
  src/quadruped_augmented.cpp
  src/quadruped_step.cpp
  src/quadruped_time.cpp
  src/quadruped_augmented_time.cpp
  src/quadruped_step_time.cpp
  )


ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE Boost::system Boost::filesystem)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC crocoddyl::crocoddyl)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include>)
IF(SUFFIX_SO_VERSION)
  SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION})
ENDIF()
INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)

ADD_SUBDIRECTORY(python)

# Build benchmark
ADD_SUBDIRECTORY(benchmark)

INSTALL(FILES package.xml DESTINATION share/${PROJECT_NAME})
