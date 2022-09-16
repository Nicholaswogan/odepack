add_library(lsoda 
  odepack_sub1.f 
  odepack_sub2.f
  odepack.f
)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(lsoda PRIVATE -std=legacy)
endif()
