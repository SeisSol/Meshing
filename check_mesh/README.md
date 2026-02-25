# Evaluate material model on mesh

## Building

```bash
mkdir build && cd build
cmake ..
make
```

## Usage

```bash
mpirun -n 4 ./check_mesh -m mesh.h5
```

## Description

1. Read `mesh.h5`
2. For each element check that: the neighbor exists
if and _only_ if the face is an internal face.
3. If an element breaks above rules, the program writes an error message.
4. If all elements are correct, the programm exits with
a success output. If at least one element is broken,
the program returns with an error code.

Internal faces are:

- Regular (0)
- Dynamic Rupture (3)
- Dynamic Rupture with Fault Tagging (64 or higher)

External faces are:

- Free Surface (1)
- Absorbing (2)
- Periodic (6)
- Free Surface with Gravity (5)
- Dirichlet (4)
- Analytical (7)
