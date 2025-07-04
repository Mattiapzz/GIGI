GIGI
=========

GIGI Generic Iterative Ggv Integration is a library specifically developed for a scientific publication submitted to RAL. This library will later evolve into a generic tool including multiple dimensions constraints and a more general approach to the problem.
The library is able to compute the maximum performance longitudinal maneuver given a generic ggv envelope and a curvature profile (i.e., the trajectory path).

## Requirements

### General

- git
- cmake
- make

### Third party libraries

Some third party libraries are needed to compile the project.

need to run the following command to install the third party libraries:

```{shell}
./third_party.sh
```

The library is self-contained and standalone, so you don't need to install any other library to use it as is. Third party libraries are used for testing and development purposes.
  
## Building

To install the project you need to run the following commands:

```{shell}
./build.sh
```

Default will build the project in release mode, if you want to build in debug mode you can run:

```{shell}
./build.sh -debug
```

For further information you can run:

```{shell}
./build.sh -h
```

### Alternative build

If you want to install the project in a different directory you can run:

```{shell}
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" ..
make -j
```
