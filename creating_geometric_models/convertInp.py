#!/usr/bin/env python3
# convert inp (from SimModeler 2d mesh (ABAQUS 2D)) to ts (Gocad)
import argparse
import numpy as np
import os
from Face import Face


def ParseInp(fname, bIsolate):
    ######## Parse input file
    with open(args.inp_filename) as fid:
        lines = fid.readlines()

    for i, line in enumerate(lines):
        if line.startswith("*Node"):
            inodes = i + 1
        if line.startswith("*Element"):
            iel = i + 1
            break
    nlines = len(lines)

    nvertex = iel - inodes - 1
    vertex = np.zeros((nvertex, 3), dtype=object)
    vid = np.zeros((nvertex), dtype=int)

    # Fill in vertex array
    vid_lookup = {}
    for i0, i in enumerate(range(inodes, iel - 1)):
        vals = lines[i].split(",")
        vid_lookup[int(vals[0])] = i0
        vertex[i0, :] = [float(val) for val in vals[1:4]]

    # Fill in the triangles array
    faces = []
    triangles = []

    for i in range(iel, nlines):
        line = lines[i]
        if line.startswith("*Element"):
            if bIsolate:
                myFace = Face(vertex=None, connect=np.asarray(triangles))
                faces.append(myFace)
                triangles = []
            continue
        val = [int(v) for v in lines[i].split(",")[1:4]]
        triangles.append(val)

    myFace = Face(vertex=None, connect=np.asarray(triangles))
    faces.append(myFace)

    for myFace in faces:
        myFace.reindex(vid_lookup)
        # reorder vertex from 0 to n where n is the number of vertice in the face
        unique_vid = list(set(list(myFace.connect.flatten())))
        vid_lu = {unique_vid[k]: k for k in range(len(unique_vid))}
        myFace.reindex(vid_lu)
        vertex0 = vertex[unique_vid, :]
        myFace.vertex = vertex0
    return faces


parser = argparse.ArgumentParser(description="convert inp (from SimModeler5 2d mesh (ABAQUS 2D)) to ts (Gocad), stl or bstl")
parser.add_argument("inp_filename", help="inp filename (SimModeler5 2d mesh (ABAQUS 2D))")
parser.add_argument("output_filename", nargs="?", help="output filname (if not used = inpbasename.ts)", default="")
parser.add_argument("--isolate", dest="isolate", action="store_true", help="isolate every surface in a different ts file")
parser.add_argument("--enforce_min_depth", nargs=1, help="move all non boundary vertices higher than zmin to zmin", type=float, metavar=("zmin"))
args = parser.parse_args()


if args.output_filename == "":
    args.output_filename = args.inp_filename[0:-4] + ".ts"
basename, ext = os.path.splitext(args.output_filename)

faces = ParseInp(args.inp_filename, args.isolate)

for i, myFace in enumerate(faces):
    if args.enforce_min_depth:
        myFace.enforce_min_depth(args.enforce_min_depth[0])

    fname = basename + str(i) + ext
    if ext in [".stl", ".ts"]:
        myFace.write(args.output_filename, append=(i != 0))
    else:
        fname = basename + str(i) + ext
        myFace.write(fname)
