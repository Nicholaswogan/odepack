cmake_minimum_required(VERSION "3.14")

project(LSODA LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

add_subdirectory(src)
add_subdirectory(test)
