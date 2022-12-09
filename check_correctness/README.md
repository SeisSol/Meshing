# Evaluate material model on mesh

## Compilation

```
mkdir build && cd build
cmake ..
make
```

## Usage

```
mpirun -n 4 ./check_correctnes -m mesh.h5
```

## Description?

1. Read `mesh.h5`
2. For each element check that:
  - if a face is a regular or a dynamic rupture face, the neighbor *has to exist*.
  - if a face is a free-surface or an absorbing face, the neighbor *must not exist*.
3. The program aborts, when the first face is found that breaks the rules above.
4. If all elements are correct, the programm exits with a succes output.
