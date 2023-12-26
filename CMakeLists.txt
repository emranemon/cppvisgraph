# The MIT License (MIT)

# Copyright (c) 2023 MIEA MD EMON 
# https://github.com/emranemon

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.10.0)

project(cppvisgraph)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find the Python interpreter and pybind11 package
## Modify the path according to the desired Python installation (If there are multiple python distributions installed)
# set(Python_EXECUTABLE "/Library/Frameworks/Python.framework/Versions/3.8/bin/python3")
find_package(Python REQUIRED COMPONENTS Interpreter Development)
## Modify the path according to the pybind11 installation (If there are multiple python distributions installed)
# set(CMAKE_PREFIX_PATH "/Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/site-packages/pybind11" ${CMAKE_PREFIX_PATH})
find_package(pybind11 REQUIRED)
# Include cppvisgraph headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
# Include Python interpreter and pybind11 headers
include_directories(${pybind11_INCLUDE_DIRS})
include_directories(${Python_INCLUDE_DIRS})
# Add the pybind11 module
pybind11_add_module(cppyvisgraph MODULE src/vg_pybind11.cpp)
# Add the cppvisgraph example module
add_executable(vgexample examples/vgexamples.cpp)
# Link against pthread library
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(vgexample PRIVATE Threads::Threads)
# Link against Python and pthread libraries
target_link_libraries(cppyvisgraph PRIVATE ${Python_LIBRARIES} Threads::Threads)