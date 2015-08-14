# Copyright (c) 2011 Daniel Pfeifer <daniel@pfeifer-mail.de>
#               2011-2014 Stefan.Eilemann@epfl.ch

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(vmmlib)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 8)
set(VERSION_PATCH 0)
set(VERSION_ABI 3)

set(COMMON_PROJECT_DOMAIN ch.eyescale)
set(DPUT_HOST "ppa:eilemann/equalizer-dev")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake
                              ${CMAKE_SOURCE_DIR}/CMake/common)
include(GitExternal)
include(Common)

common_package(Boost COMPONENTS unit_test_framework)
common_package_post()

set(DOCS README.md AUTHORS LICENSE.txt ACKNOWLEDGEMENTS)

add_subdirectory(vmmlib)
add_subdirectory(tests)
add_subdirectory(tests_old)

install(FILES ${DOCS} DESTINATION share/vmmlib COMPONENT dev)

include(DoxygenRule) # must be after all targets
include(CPackConfig)
include(CTest)
