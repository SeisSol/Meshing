#!/usr/bin/env python3
from Face import Face
import argparse
import os
import trimesh
import numpy as np

parser = argparse.ArgumentParser(description="refine and smooth")
parser.add_argument("input_file", help="surface in ts file format")
parser.add_argument(
    "--N",
    type=int,
    required=True,
    help="number of initial refine steps",
)
parser.add_argument(
    "--P",
    type=int,
    required=True,
    help="number of and refine/smoothing steps",
)

parser.add_argument(
    "--remesh_first",
    nargs=1,
    metavar=("mesh_size"),
    type=float,
    help="remesh first the surface with pygalmesh",
)

parser.add_argument(
    "--fix_boundary",
    dest="fix_boundary",
    action="store_true",
    help="fix boundary edge when interpolating",
)


args = parser.parse_args()

myFace = Face.from_file(args.input_file)

if args.remesh_first:
    """
    # pygalmesh installed with:
    sudo apt install libcgal-dev libeigen3-dev
    pip install pygalmesh
    """
    import pygalmesh
    import meshio

    mesh = meshio.Mesh(points=myFace.vertex, cells=[("triangle", myFace.connect)])
    mesh.write("trash_me.vtk")
    mesh_size = args.remesh_first[0]
    print(f"remeshing with pygalmesh aiming for mesh size {mesh_size}")
    mesh = pygalmesh.remesh_surface(
        "trash_me.vtk",
        max_edge_size_at_feature_edges=mesh_size,
        min_facet_angle=25,
        max_radius_surface_delaunay_ball=mesh_size,
        max_facet_distance=mesh_size,
        verbose=False,
    )
    os.remove("trash_me.vtk")
    print("done remeshing")
    a = trimesh.Trimesh(vertices=mesh.points, faces=mesh.cells[0].data)
else:
    a = trimesh.Trimesh(vertices=myFace.vertex, faces=myFace.connect)

for i in range(args.N):
    a = a.subdivide()
for i in range(args.P):
    a = a.subdivide()

    if args.fix_boundary:
        # Identify boundary edges using grouping
        unique_edges = a.edges[
            trimesh.grouping.group_rows(a.edges_sorted, require_count=1)
        ]
        # Extract unique vertices on boundary
        boundary_vertices = np.unique(unique_edges)
        a_before_smoothing = a.copy()

    a = trimesh.smoothing.filter_taubin(a)

    if args.fix_boundary:
        # Replace only the internal vertices
        a.vertices[np.isin(np.arange(len(a.vertices)), boundary_vertices)] = (
            a_before_smoothing.vertices[
                np.isin(np.arange(len(a.vertices)), boundary_vertices)
            ]
        )


myFace = Face(a.vertices, a.faces)
basename, ext = os.path.splitext(args.input_file)
myFace.write(f"{basename}_refined_smooth_{args.N}.ts")
