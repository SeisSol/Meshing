gmsh2gambit
===============
Transform msh file (gmsh msh 2.0 mesh file) to neu file (neutral mesh file)
From gmsh 4.0, the option '-format neu' of gmsh should be used instead.

Build
===============
Using CMake (preferred)
```
cmake
make -j 4
```

Using scons (legacy)
```
CC=gcc CXX=g++ scons   
```
Or, replace `gcc` and `g++` by another compiler (e.g. `icc` and `icpc` for the legacy Intel compilers).
Note that the specification of `CC` and `CXX` is mandatory for the scons script to work.

Usage
===============
gmsh2gambit -i test.msh -o test.neu


