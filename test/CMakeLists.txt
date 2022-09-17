set(ODEPACK_TESTS test_lsoda test_threadsafe benchmark)

foreach(ODEPACK_TEST ${ODEPACK_TESTS})
  add_executable(${ODEPACK_TEST} "${ODEPACK_TEST}.f90")
  target_link_libraries(${ODEPACK_TEST} odepack_mod)
  target_include_directories(${ODEPACK_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
endforeach()

# multithreading flags
find_package(OpenMP REQUIRED)
if (OpenMP_Fortran_FOUND)
  target_compile_options(test_threadsafe PRIVATE ${OpenMP_Fortran_FLAGS})
  target_link_options(test_threadsafe PUBLIC ${OpenMP_Fortran_FLAGS})
endif()