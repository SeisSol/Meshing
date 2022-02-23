#!/usr/bin/env python3
import pygmsh
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="generate a meshed plane with pygmsh")
parser.add_argument("output_file", help="filename of output file")

parser.add_argument(
    "--center",
    nargs=3,
    metavar=("x0", "y0", "z0"),
    help="center of rectangle",
    type=float,
    required=True,
)

parser.add_argument(
    "--normal",
    nargs=3,
    metavar=("nx", "ny", "nz"),
    help="normal orientation",
    type=float,
    required=True,
)

parser.add_argument(
    "--dims",
    nargs=2,
    metavar=("lx", "lz"),
    help="plane dimensions",
    type=float,
    required=True,
)

parser.add_argument(
    "--meshSize",
    nargs=1,
    metavar=("size"),
    help="2d mesh size",
    type=float,
    default=[10e3],
)

args = parser.parse_args()

nx, ny, nz = args.normal
normal = np.array(args.normal)
if max(abs(nx), abs(ny)) == 0:
    ux = np.array([1.0, 0, 0])
    uz = np.array([0, 1.0, 0])
ux = np.array([-ny, nx, 0.0])
uz = np.cross(normal, ux)
ux = ux / np.linalg.norm(ux)
uz = uz / np.linalg.norm(uz)

center = np.array(args.center)
lx, lz = args.dims
p0 = center - 0.5 * lx * ux - 0.5 * lz * uz
p1 = center - 0.5 * lx * ux + 0.5 * lz * uz
p2 = center + 0.5 * lx * ux + 0.5 * lz * uz
p3 = center + 0.5 * lx * ux - 0.5 * lz * uz

with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_max = args.meshSize[0]
    geom.add_polygon(
        [p0, p1, p2, p3],
    )
    mesh = geom.generate_mesh(dim=2)
    mesh.write(args.output_file)
