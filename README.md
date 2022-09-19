# ODEPACK

This is a Modern Fortran interface the `LSODA` and `LSODAR` routines in [ODEPACK](https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html), which is for solving ordinary differential equation initial value problems. This repository contains a modified version of ODEPACK which is threadsafe.

## Example

For a simple example see `test/test_lsoda.f90`

## Building

**CMake**

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
# run tests
./test/test_lsoda
./test/test_threadsafe
# run benchmark
./test/benchmark
```

Add as a dependency to your CMake projects with the [CMake Package Manager](https://github.com/cpm-cmake/CPM.cmake)

**Fortran Package Manager**

```sh
fpm build --profile release
fpm test --profile release
```

## Documentation

The file `src/odepack_mod.f90` contains extensive comments describing the use of the user-interfacing class `lsoda_class`.

## Differences from original ODEPACK

- This repository replaces `common` blocks with a passed derived type. This makes the code threadsafe.
- This repository changes the interface to the ODE right-hand-side subroutine, the jacobian and the root-finding subroutine, by adding a `ierr` argument. If `ierr` is set to less than 0 in these subroutines, then the integration is gracefully terminated.
- Root-finding now indicates the direction of the root. See variable `jroot` in `src/odepack_mod.f90`