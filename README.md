# ODEPACK

This is a Modern Fortran interface the `LSODA` and `LSODAR` routines in [ODEPACK](https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html). This repository contains a modified version of ODEPACK which is threadsafe.

## Example

For a simple example see `test/test_lsoda.f90`

## Building

### CMake

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

### Fortran Package Manager

```sh
fpm build --profile release
fpm test --profile release
```