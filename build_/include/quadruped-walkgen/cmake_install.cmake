# Install script for directory: /home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/tcorberes/Desktop/sfe/quadruped_walkgen_python")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/quadruped-walkgen" TYPE FILE FILES
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_nl.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_nl.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_augmented.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_augmented.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_step.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_step.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_augmented_time.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_augmented_time.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_step_time.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_step_time.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_step_period.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_step_period.hxx"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_time.hpp"
    "/home/tcorberes/Desktop/sfe/quadruped-walkgen/include/quadruped-walkgen/quadruped_time.hxx"
    )
endif()

