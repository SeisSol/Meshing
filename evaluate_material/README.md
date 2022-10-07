# Evaluate material model on mesh

## Compilation

```
mkdir build && cd build
cmake ..
make
```

## Usage

```
./evaluate_easi -m mesh.h5 -e easi.yaml -o material
```

## Description?

1. Compute barycenter for each tetrahedron in `mesh.h5`.
2. Evaluate all parameters, which are described in `easi.yaml`,b at the element barycenters.
3. Write `material.xdmf` and `material.h5`, which can be visualized with paraview.
