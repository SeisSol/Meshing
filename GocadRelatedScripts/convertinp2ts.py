# convert inp (from SimModeler 2d mesh (ABAQUS 2D)) to ts (Gocad)

import argparse
import numpy as np
from collections import defaultdict


class Face:
    def __init__(self, connect):
        self.connect = connect
        self.id_vertex = np.unique(self.connect.flatten())
        self.local_vid_lookup = {}
        for i, ivx in enumerate(self.id_vertex):
            self.local_vid_lookup[ivx] = i
        self.ntriangles = self.connect.shape[0]

    def intersect(self, Face2):
        "return the common nodes between self and Face2"
        return np.intersect1d(self.id_vertex, Face2.id_vertex)

    def writeTs(self, fname, vid_lookup, vertex, control_nodes=[]):
        "output face as a *.ts file"
        iscontroleNode = np.zeros((len(myFace.id_vertex)), dtype=int)
        for cn in control_nodes:
            iscontroleNode[self.local_vid_lookup[cn]] = 1
        with open(fname, "w") as fout:
            fout.write("GOCAD TSURF 1\nHEADER {\nname:%s\nborder: true\nmesh: false\n*border*bstone: true\n}\nTFACE\n" % (fname))
            for kk, ivx in enumerate(myFace.id_vertex):
                i = vid_lookup[ivx]
                if iscontroleNode[kk]:
                    fout.write("VRTX %s %s %s %s CNXYZ\n" % (ivx, vertex[i, 0], vertex[i, 1], vertex[i, 2].rstrip()))
                else:
                    fout.write("VRTX %s %s %s %s" % (ivx, vertex[i, 0], vertex[i, 1], vertex[i, 2]))

            for i in range(self.ntriangles):
                fout.write("TRGL %d %d %d\n" % (self.connect[i, 0], self.connect[i, 1], self.connect[i, 2]))

            fout.write("END\n")

        print("done writing " + fname)


parser = argparse.ArgumentParser(description="convert inp (from SimModeler5 2d mesh (ABAQUS 2D)) to ts (Gocad)")
parser.add_argument("inp_filename", help="inp filename (SimModeler5 2d mesh (ABAQUS 2D))")
parser.add_argument("ts_filename", nargs="?", help="output filname (if not used = inpbasename.ts)", default="")
parser.add_argument("--isolate", dest="isolate", action="store_true", help="isolate every surface in a different ts file")
parser.add_argument("--tag_intersect", dest="tag_intersect", action="store_true", help="write the intersection in the ts file as control nodes")
args = parser.parse_args()

if args.ts_filename == "":
    args.ts_filename = args.inp_filename[0:-4] + ".ts"

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
for i in range(inodes, iel - 1):
    vals = lines[i].split(",")
    vid_lookup[int(vals[0])] = i - inodes
    vertex[i - inodes, :] = vals[1:4]

# Fill in the triangles array
faces = []
triangles = []

for i in range(iel, len(lines)):
    line = lines[i]
    if line.startswith("*Element"):
        if args.isolate:
            triangles = np.asarray(triangles)
            myFace = Face(triangles)
            faces.append(myFace)
            triangles = []
        continue
    val = [int(v) for v in lines[i].split(",")[1:4]]
    triangles.append(val)

myFace = Face(np.asarray(triangles))
faces.append(myFace)

if args.tag_intersect:
    print("generating intersection...")
    from itertools import combinations

    inters = {}
    for t1, t2 in combinations(range(len(faces)), 2):
        inters[tuple((t1, t2))] = faces[t1].intersect(faces[t2])

for i, myFace in enumerate(faces):
    fname = args.ts_filename[0:-3] + str(i) + ".ts"

    if args.tag_intersect:
        myinters = []
        for t12 in inters.keys():
            if i in t12:
                myinters.extend(inters[t12])
        myFace.writeTs(fname, vid_lookup, vertex, myinters)
    else:
        myFace.writeTs(fname, vid_lookup, vertex)
