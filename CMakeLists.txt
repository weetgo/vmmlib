# Copyright (c) 2011-2016 Daniel Pfeifer   daniel@pfeifer-mail.de
#                         Steafan Eilemann Stefan.Eilemann@epfl.ch

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(vmmlib VERSION 1.12.0)
set(vmmlib_VERSION_ABI 7)

set(COMMON_PROJECT_DOMAIN ch.eyescale)
set(VMMLIB_DESCRIPTION "Templatized C++ vector and matrix math library")
set(VMMLIB_MAINTAINER "Stefan Eilemann <eile@eyescale.ch>")
set(VMMLIB_LICENSE BSD)
set(DPUT_HOST "ppa:eilemann/equalizer-dev")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
                              ${CMAKE_SOURCE_DIR}/CMake/common)
include(GitExternal)
include(Common)

common_find_package(Boost COMPONENTS unit_test_framework)
common_find_package_post()

add_subdirectory(vmmlib)
add_subdirectory(tests)

set(DOCS README.md LICENSE.txt ACKNOWLEDGEMENTS)
install(FILES ${DOCS} DESTINATION share/vmmlib COMPONENT dev)

set(DOXYGEN_MAINPAGE_MD README.md)
set(DOXYGEN_EXTRA_INPUT ${PROJECT_SOURCE_DIR}/README.md ${PROJECT_SOURCE_DIR}/doc)
include(DoxygenRule) # must be after all targets
include(CPackConfig)
include(CTest)
