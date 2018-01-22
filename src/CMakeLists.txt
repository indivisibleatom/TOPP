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

# copy the python files
set(CMAKE_BUILD_TYPE Debug)

# include python bindings sample
if( OpenRAVE_FOUND )
  include_directories(${PYTHON_INCLUDE_PATH} ${OpenRAVE_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
  link_directories(${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})
else()
  include_directories(${PYTHON_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
endif()
file(GLOB topp_files ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(TOPPbindings SHARED ${topp_files})
add_python_module(TOPP python TOPP)

# Output bindings to python module dir
set_target_properties(TOPPbindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/python/TOPP)

if( OpenRAVE_FOUND )
  target_link_libraries(TOPPbindings ${OpenRAVE_PYTHON_DIR}/openravepy_int.so ${PYTHON_LIBS} ${Boost_THREAD_LIBRARY} ${OpenRAVE_LIBRARIES} ${LAPACK_LIBRARIES})
  set_target_properties(TOPPbindings PROPERTIES PREFIX "" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} ${OpenRAVE_CXX_FLAGS}" LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")
else()
  target_link_libraries(TOPPbindings ${PYTHON_LIBS} ${Boost_THREAD_LIBRARY} ${LAPACK_LIBRARIES})
  set_target_properties(TOPPbindings PROPERTIES PREFIX "" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} " LINK_FLAGS "")
endif()

