# Copyright (c) 2011-2016 Daniel Pfeifer <daniel@pfeifer-mail.de>
#                         Stefan.Eilemann@epfl.ch

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(vmmlib)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 10)
set(VERSION_PATCH 0)
set(VERSION_ABI 5)

set(COMMON_PROJECT_DOMAIN ch.eyescale)
set(VMMLIB_DESCRIPTION "Templatized C++ vector and matrix math library")
set(VMMLIB_MAINTAINER "Stefan Eilemann <eile@eyescale.ch>")
set(VMMLIB_LICENSE BSD)
set(DPUT_HOST "ppa:eilemann/equalizer-dev")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
                              ${CMAKE_SOURCE_DIR}/CMake/common)
include(GitExternal)
include(Common)

common_package(Boost COMPONENTS unit_test_framework)
common_package_post()

add_subdirectory(vmmlib)
add_subdirectory(tests)
add_subdirectory(tests_old)

set(DOCS README.md LICENSE.txt ACKNOWLEDGEMENTS)
install(FILES ${DOCS} DESTINATION share/vmmlib COMPONENT dev)

set(DOXYGEN_MAINPAGE_MD README.md)
set(DOXYGEN_EXTRA_INPUT ${PROJECT_SOURCE_DIR}/README.md)
include(DoxygenRule) # must be after all targets
include(CPackConfig)
include(CTest)
