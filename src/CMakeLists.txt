add_library(odepack
  odepack_interface.f90
  odepack_common.f90
  odepack_sub1.f 
  odepack_sub2.f
  odepack.f
)

if (${REQUIRE_FAST_LAPACK})
  find_package(LAPACK REQUIRED)
else()
  find_package(LAPACK)
endif()
if (${LAPACK_FOUND})
  target_link_libraries(odepack ${LAPACK_LIBRARIES})
else()
  # We compile lapack files downloaded from netlib
  file(GLOB LAPACK_SRC_FILES lapack/*.f)
  add_library(lapack ${LAPACK_SRC_FILES})
  target_link_libraries(odepack lapack)
endif()

add_library(odepack_mod
  odepack_mod.f90
)
target_link_libraries(odepack_mod odepack)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(odepack PRIVATE -std=legacy)
  if ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    target_compile_options(odepack PRIVATE -funroll-loops)
    target_compile_options(odepack_mod PRIVATE -funroll-loops)
  endif()
endif()
