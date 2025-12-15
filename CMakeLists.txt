cmake_minimum_required(VERSION "3.14")

project(LSODA LANGUAGES Fortran VERSION "0.2.1")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

option(ENABLE_PYTHON_TESTS "Build Python test driver targets" OFF)

option(REQUIRE_FAST_LAPACK 
"If ON, then cmake configuration will 
error out if e.g. OpenBLAS is not found" OFF)

add_subdirectory(src)
add_subdirectory(test)
if (ENABLE_PYTHON_TESTS)
  add_subdirectory(python/tests)
endif()
