Place receivers
===============
Given receiver locations on a 2D plane, this tool finds the corresponding elevations on the free-surface of the given mesh.

Method
------
The program iterates over all faces with free surface boundary condition. For each face, the tool checks if a receiver lies inside the projection of the face on the x-y-plane. Then, the receiver's z coordinate is set such that the receiver lies on the face. Afterwards it is lowered by a configurable amount along the z axis. (The direction is determined by the face normal, i.e. "lowered" means in the opposite direction of the normal.)

Compilation 
-----

```
mkdir build && cd build
cmake ..
make
```

Usage
-----

```
place_receivers -d DEPTH -r RECEIVERS -m MESH -o OUTPUT
```

- DEPTH should be a positive number that specifies the amount of meters that the receiver shall be placed below the free surface.
- RECEIVERS is the file name of a double column, space separated file where the first column denotes the x coordinate and the second column the y coordinate. For example:

<pre>
8.251411e+05 5.900044e+05
7.654161e+05 3.274531e+05
7.794138e+05 5.097186e+05
</pre>

One may specify the depth per receiver with an optional third column, where the depth is the maximum of DEPTH and the third column. For example with DEPTH=0.1:

<pre>
8.251411e+05 5.900044e+05 0.0   # -> depth 0.1
7.654161e+05 3.274531e+05 100.0 # -> depth 100
</pre>

- MESH is the file name of the mesh.
- OUTPUT is the file name of the output file.

Limitations
-----------
The tools does only work if faces with free surface boundary condition are present and if their normals look in +z or -z direction. (This is usually the case if free surfaces due to topography and e.g. enu, ned, seu axis orientations are used.) That is, the tool fails if all respective faces have normals with z=0.
