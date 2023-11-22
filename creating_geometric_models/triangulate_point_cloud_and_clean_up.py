#!/usr/bin/env python3
from sklearn.decomposition import PCA
import numpy as np
from scipy.spatial import Delaunay
import trimesh
from Face import Face
import argparse
import os

parser = argparse.ArgumentParser(
    description=(
        "triangulate a point cloud and recursively clean the triangulation based on"
        " edge-length, area and aspect ratio criteria"
    )
)
parser.add_argument("input_file", help="input file (3 columns ascii x y z)")
parser.add_argument(
    "--edge_length_threshold",
    nargs=1,
    type=float,
    help="border triangles with edge length larger than the threshold will be removed",
)
parser.add_argument(
    "--area_threshold",
    nargs=1,
    type=float,
    help="border triangles with area lower than the threshold will be removed",
)
parser.add_argument(
    "--aspect_ratio_threshold",
    nargs=1,
    type=float,
    default=[7.0],
    help="border triangles with aspect ratio larger than the threshold will be removed",
)

args = parser.parse_args()


xyz = np.loadtxt(args.input_file)
# Perform PCA to get principal axes
pca = PCA(n_components=2)
points = pca.fit_transform(xyz)
tri = Delaunay(points)


mesh = trimesh.Trimesh(vertices=xyz, faces=tri.simplices)
edge_lengths = mesh.edges_unique_length

edge_length_threshold = (
    np.percentile(edge_lengths, 96)
    if not args.edge_length_threshold
    else args.edge_length_threshold[0]
)
print(f"using edge_length_threshold = {edge_length_threshold}")

triangle_areas = mesh.area_faces
triangle_area_threshold = (
    np.percentile(triangle_areas, 3)
    if not args.area_threshold
    else args.area_threshold[0]
)
print(f"using triangle_area_threshold = {triangle_area_threshold}")

print(f"using aspect_ratio_threshold = {args.aspect_ratio_threshold[0]}")

k = 0
while True:
    index = trimesh.grouping.group_rows(mesh.edges_sorted, require_count=1)
    boundary_faces_id = mesh.edges_face[index]
    triangle_areas = mesh.area_faces[boundary_faces_id]

    edge_lengths = np.zeros(len(boundary_faces_id))
    aspect_ratio = np.zeros(len(boundary_faces_id))
    for i, face in enumerate(mesh.faces[boundary_faces_id]):
        edge_lengths3 = np.linalg.norm(
            mesh.vertices[face] - np.roll(mesh.vertices[face], shift=-1, axis=0), axis=1
        )
        edge_lengths[i] = np.max(edge_lengths3)
        aspect_ratio[i] = np.max(edge_lengths3) / np.min(edge_lengths3)

    ar_condition = aspect_ratio > args.aspect_ratio_threshold[0]
    length_condition = edge_lengths > edge_length_threshold
    area_condition = triangle_areas < triangle_area_threshold

    combined_condition = np.logical_or(ar_condition, length_condition)
    combined_condition = np.logical_or(combined_condition, area_condition)
    index = index[combined_condition]

    index_faces = mesh.edges_face[index]
    if len(index_faces) == 0:
        break
    mesh = mesh.submesh([~np.isin(np.arange(len(mesh.faces)), index_faces)])[0]
    print(
        f"iteration {k}, {len(index_faces)} triangles removed,"
        f" {len(mesh.faces)} remaining"
    )
    k = k + 1

myFace = Face(mesh.vertices, mesh.faces)
basename, ext = os.path.splitext(args.input_file)
myFace.write(f"{basename}_delaunay_clean.ts")
