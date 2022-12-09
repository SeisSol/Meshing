# Evaluate material model on mesh

## Compilation

```
mkdir build && cd build
cmake ..
make
```

## Usage

```
mpirun -n 4 ./check_mesh -m mesh.h5
```

## Description

1. Read `mesh.h5`
2. For each element check that:
  - if a face is an internal face, the neighbor *has to exist*.
  - if a face is an external face, the neighbor *must not exist*.
3. If an element breaks above rules, the program writes an error message.
3. If all elements are orrect, the programm exits with a success output. If at least one element is broken, the program returns with an error code.

Internal faces are:
* Regular
* Dynamic Rupture
* Dynamic Rupture with Fault Tagging
External Faces are:
* Free Surface
* Absorbing
* Periodic
* Free Surface with Gravity
