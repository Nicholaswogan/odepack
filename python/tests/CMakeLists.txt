# Fortran driver used by Python tests (writes binary reference outputs).
add_executable(odepack_sub1_driver "${CMAKE_CURRENT_SOURCE_DIR}/test_odepack_sub1.f90")
target_link_libraries(odepack_sub1_driver odepack_mod)
target_include_directories(odepack_sub1_driver PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(odepack_sub1_driver PRIVATE -std=legacy)
endif()
