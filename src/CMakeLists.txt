add_library(odepack
  odepack_interface.f90
  odepack_common.f90
  odepack_sub1.f 
  odepack_sub2.f
  odepack.f
)

add_library(odepack_mod
  odepack_mod.f90
)
target_link_libraries(odepack_mod odepack)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(odepack PRIVATE -std=legacy -w)
endif()
