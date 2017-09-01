# -*- coding: utf-8 -*-
# Copyright (C) 2014 Rosen Diankov <rosen.diankov@gmail.com> & Cuong Pham
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
cmake_minimum_required (VERSION 2.8.0)
project(topp)
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

set (TOPP_VERSION_MAJOR 0)
set (TOPP_VERSION_MINOR 1)
set (TOPP_VERSION_PATCH 0)
set (TOPP_VERSION ${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR}.${TOPP_VERSION_PATCH})
set (TOPP_SOVERSION ${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR})
set (CLIENT_SOVERSION ${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR})
message(STATUS "Compiling TOPP Version ${TOPP_VERSION}, soversion=${CLIENT_SOVERSION}")

set(TOPP_CMAKE_INSTALL_DIR "topp-${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR}" CACHE STRING "Directory to install the cmake config files.")
set(TOPP_TARGET_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR} CACHE STRING "The target processor architecture to build for, this is combined with the generator toolchain")

message(STATUS "Using cmake version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" )
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0002
cmake_policy(SET CMP0002 NEW)
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0003
cmake_policy(SET CMP0003 NEW)

# Use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# When building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(LIB_SUFFIX CACHE STRING "suffix for the library directory need for x86-64 systems that use lib64 ")

# The RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(CheckCXXCompilerFlag)

# have to include before boost since the boost headers can be located in a previous installed version of this library
set(TOPP_INCLUDE_LOCAL_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${TOPP_INCLUDE_LOCAL_DIRS})

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG )
  set(CMAKE_CXX_FLAGS_OPTIMIZED "-O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS -D_SECURE_SCL=0") # this practically removes all checks making it a very dangerous options to play with
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -D_DEBUG")

  add_definitions(-fPIC) # this is a very important switch and some libraries seem now to have it....
  # for some reason c++-0x is required...
  set(EXTRA_COMPILE_FLAGS "-fPIC -std=c++0x")
else()
  set(EXTRA_COMPILE_FLAGS "")
endif()

if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )


find_package(OpenRAVE)
if( OpenRAVE_FOUND )
  set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DWITH_OPENRAVE")
endif()

#find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams python thread filesystem system thread date_time)# regex)
find_package(Boost COMPONENTS iostreams python3 thread filesystem system thread date_time)
if( Boost_FOUND )
  include_directories(${Boost_INCLUDE_DIRS})
  set(TOPP_LINK_DIRS ${TOPP_LINK_DIRS} ${Boost_LIBRARY_DIRS})
elseif(Boost_VERSION AND NOT "${Boost_VERSION}" STREQUAL "0")
  include_directories(${Boost_INCLUDE_DIRS})
  set(TOPP_LINK_DIRS ${TOPP_LINK_DIRS} ${Boost_LIBRARY_DIRS})
  message(FATAL_ERROR "Could not find boost libraries!")
endif()

message(STATUS "found boost version: ${Boost_VERSION} at ${Boost_INCLUDE_DIRS}")
if( NOT Boost_THREAD_LIBRARY OR NOT Boost_DATE_TIME_LIBRARY )
  message(FATAL_ERROR "boost thread and date_time libraries are necessary")
endif()

check_library_exists(lapack cgesdd_ "" LAPACK_FOUND)
if( LAPACK_FOUND )
  set(LAPACK_LIBRARIES lapack)
else()
  message(STATUS "lapack dev library not found")
endif()

add_subdirectory(src)
