#!/usr/bin/env python3
import numpy as np
import argparse
import os
from Face import Face
from Grid import Grid

parser = argparse.ArgumentParser(description="create surface from a (possibly sparse, e.g. Slab2.0 dataset) structured dataset (e.g. netcdf)")
parser.add_argument("input_file", help="netcdf file")
parser.add_argument("output_file", help="output file (ext in stl, bstl, ts)")
parser.add_argument("--subsample", nargs=1, type=int, metavar=("onesample_every"), default=[1], help="use only one value every onesample_every in both direction")
parser.add_argument("--objectname", nargs=1, metavar=("objectname"), default=(""), help="name of the surface in gocad")
parser.add_argument("--hole", nargs=4, metavar=(("x0"), ("x1"), ("y0"), ("y1")), help="isolate a hole in surface defined by x0<=x<=x1 and y0<=y<=y1 (stl and ts output only)", type=float)
parser.add_argument("--crop", nargs=4, metavar=(("x0"), ("x1"), ("y0"), ("y1")), help="select only surfaces in x0<=x<=x1 and y0<=y<=y1", type=float)
parser.add_argument("--proj", nargs=1, metavar=("projname"), help="transform vertex array to projected system.\
 projname: name of the (projected) Coordinate Reference System (CRS) (e.g. EPSG:32646 for UTM46N)")
parser.add_argument("--translate", nargs=2, metavar=("x0", "y0"), default=([0, 0]), help="translates all nodes by (x0,y0)", type=float)
args = parser.parse_args()

if args.objectname == "":
    base = os.path.basename(args.input_file)
    args.objectname = os.path.splitext(base)[0]
else:
    args.objectname = args.objectname[0]

structured_grid = Grid(args.input_file, args.subsample[0])
structured_grid.crop(args.crop)
structured_grid.generate_vertex()

structured_grid.generate_connect()
structured_grid.isolate_hole(args.hole)
if args.proj:
    structured_grid.proj_vertex(args.proj[0])

if args.translate:
    structured_grid.translate(args.translate)

basename, ext = os.path.splitext(args.output_file)
nsolid = max(structured_grid.solid_id) + 1


if nsolid == 1:
    myFace = Face(structured_grid.vertex, structured_grid.connect)
    if structured_grid.is_sparse:
        myFace.reindex(structured_grid.vid_lookup)
    myFace.write(f"{basename}{ext}")
else:
    for sid in range(nsolid):
        idtr = np.where(structured_grid.solid_id == sid)[0]
        aVid = np.unique(structured_grid.connect[idtr, :].flatten())
        myFace = Face(vertex=None, connect=structured_grid.connect[idtr, :])
        if structured_grid.is_sparse:
            myFace.reindex(structured_grid.vid_lookup)
        myFace.write(f"{basename}{sid}{ext}", structured_grid.vertex, write_full_vertex_array=False)
