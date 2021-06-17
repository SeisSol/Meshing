import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="create surface from an partially structured grid of nodes.\
 the point set should consist of serveral lines of nodes of same y (or x, depending on args.axis) coordinates")
parser.add_argument("input_file", help="x y z, one node coordinate by line")
parser.add_argument("output_file", help="gocad ts output file")
parser.add_argument("--axis", nargs=1, metavar=("axis"), default=(0), help="axis of structured coordinate (0:X, 1:Y)", type=int, choices=[0,1])
parser.add_argument("--objectname", nargs=1, metavar=("objectname"), help="name of the surface in gocad")
parser.add_argument("--large_precision", dest="large_precision", action="store_true", help="write vertex using \%25.20f")
parser.add_argument("--append", dest="append", action="store_true", help="append to output_file")
args = parser.parse_args()


if not args.objectname:
    base = os.path.basename(args.input_file)
    args.objectname = os.path.splitext(base)[0]
else:
    args.objectname = args.objectname[0]

inputfn = args.input_file
axis = args.axis[0]

dataxyz = np.loadtxt(inputfn)
nvertex = np.shape(dataxyz)[0]

structured_coord = set(dataxyz[:, axis])
structured_coord = sorted(structured_coord)
# precision for selection the nodes with a given latitude
dx = 1e-4
triangles = []

for i, val in enumerate(structured_coord):
    if i == 0:
        p1 = np.where(abs(dataxyz[:, axis] - val) < dx)[0]
        np1 = np.size(p1)
        continue
    p0 = p1
    np0 = np1
    p1 = np.where(abs(dataxyz[:, axis] - val) < dx)[0]
    np1 = np.size(p1)
    if np0 < np1:
        for i0 in range(np0 - 1):
            triangles.append([p0[i0], p0[i0 + 1], p1[i0]])
            triangles.append([p0[i0 + 1], p1[i0 + 1], p1[i0]])
        for i0 in range(np0 - 1, np1 - 1):
            triangles.append([p0[np0 - 1], p1[i0 + 1], p1[i0]])
    else:
        for i1 in range(np1 - 1):
            triangles.append([p1[i1], p0[i1], p1[i1 + 1]])
            triangles.append([p1[i1 + 1], p0[i1], p0[i1 + 1]])
        for i1 in range(np1 - 1, np0 - 1):
            triangles.append([p1[np1 - 1], p0[i1], p0[i1 + 1]])

if args.append:
    fout = open(args.output_file, "a")
else:
    fout = open(args.output_file, "w")
fout.write("GOCAD TSURF 1\nHEADER {\nname:" + args.objectname + "\n}\nTRIANGLES\n")

for i in range(0, nvertex):
    if args.large_precision:
        fout.write("VRTX %d %25.20f %25.20f %25.20f\n" % (i + 1, dataxyz[i, 0], dataxyz[i, 1], dataxyz[i, 2]))
    else:
        fout.write("VRTX %d %f %f %f\n" % (i + 1, dataxyz[i, 0], dataxyz[i, 1], dataxyz[i, 2]))

for tr in triangles:
    fout.write("TRGL %d %d %d\n" % (tr[0] + 1, tr[1] + 1, tr[2] + 1))


fout.write("END\n")
