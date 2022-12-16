cmake_minimum_required(VERSION 3.14)

# Set up project properties
set(PROJECT_NAMESPACE loco-3d)
set(PROJECT_NAME quadruped-walkgen)
set(PROJECT_DESCRIPTION "MPC Controller using Crocoddyl")
set(PROJECT_URL "https://github.com/${PROJECT_NAMESPACE}/${PROJECT_NAME}")

option(SUFFIX_SO_VERSION "Suffix library name with its version" ON)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CUSTOM_HEADER_DIR ${PROJECT_NAME})
set(CXX_DISABLE_WERROR TRUE)
set(PROJECT_USE_CMAKE_EXPORT TRUE)

# Check if the submodule cmake have been initialized Check if the submodule
# cmake have been initialized
set(JRL_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/cmake")
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/cmake/base.cmake")
  message(STATUS "JRL cmakemodules not found. Let's fetch it.")
  include(FetchContent)
  FetchContent_Declare(
    "jrl-cmakemodules"
    GIT_REPOSITORY "https://github.com/jrl-umi3218/jrl-cmakemodules.git")
  FetchContent_MakeAvailable("jrl-cmakemodules")
  FetchContent_GetProperties("jrl-cmakemodules" SOURCE_DIR JRL_CMAKE_MODULES)
endif()

# Include important cmake modules
include("${JRL_CMAKE_MODULES}/base.cmake")

compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})
# ADD_DEFINITIONS(-DWITH_CODEGEN)

add_project_dependency(crocoddyl REQUIRED)
add_project_dependency(example-robot-data)
find_package(Boost REQUIRED COMPONENTS filesystem system)
string(REGEX REPLACE "-" "_" PYTHON_DIR ${PROJECT_NAME})

option(BUILD_WITH_MULTITHREADS
       "Build the library with the Multithreading support (required OpenMP)" ON)
if(BUILD_WITH_MULTITHREADS)
  set(BUILD_WITH_NTHREADS
      "4"
      CACHE STRING "Number of threads")
  string(REGEX MATCH "^[0-9]+$" BUILD_WITH_NTHREADS ${BUILD_WITH_NTHREADS})
  if(NOT BUILD_WITH_NTHREADS MATCHES "^[0-9]+$")
    set(BUILD_WITH_NTHREADS 4)
    message(
      WARNING
        "the number of threads have to be an interger value, set to ${BUILD_WITH_NTHREADS}"
    )
  endif()
endif()

set(${PROJECT_NAME}_HEADERS
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
    include/${CUSTOM_HEADER_DIR}/quadruped_time.hxx)

set(${PROJECT_NAME}_SOURCES
    src/quadruped.cpp
    src/quadruped_nl.cpp
    src/quadruped_augmented.cpp
    src/quadruped_step.cpp
    src/quadruped_time.cpp
    src/quadruped_augmented_time.cpp
    src/quadruped_step_time.cpp)

add_library(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES}
                                   ${${PROJECT_NAME}_HEADERS})
target_link_libraries(${PROJECT_NAME} PRIVATE Boost::system Boost::filesystem)
target_link_libraries(${PROJECT_NAME} PUBLIC crocoddyl::crocoddyl)
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include>)
if(SUFFIX_SO_VERSION)
  set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION})
endif()
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${TARGETS_EXPORT_NAME}
  DESTINATION lib)

add_subdirectory(python)

# Build benchmark
add_subdirectory(benchmark)

install(FILES package.xml DESTINATION share/${PROJECT_NAME})
