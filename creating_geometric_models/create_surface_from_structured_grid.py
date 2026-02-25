#!/usr/bin/env python3
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="create surface from a structured grid of nodes.\
 vertices in a row or a column do not necessary share the same x and y values (non rectilinear grid)")
parser.add_argument("input_file", help="x y z, one node coordinate by line")
parser.add_argument("output_file", help="gocad output file")
parser.add_argument("--NX", nargs=1, metavar=("NX"), help="NX: number of nodes in the first structured dimension")
parser.add_argument("--subsample", nargs=1, metavar=("onesample_every"), help="use only one value every onesample_every in both direction")
parser.add_argument("--objectname", nargs=1, metavar=("objectname"), help="name of the surface in gocad")
parser.add_argument("--hole", nargs=4, metavar=(("x0"), ("x1"), ("y0"), ("y1")), help="create a hole in surface defined by x0<=x<=x1 and y0<=y<=y1")
parser.add_argument("--crop", nargs=4, metavar=(("x0"), ("x1"), ("y0"), ("y1")), help="select only surfaces in x0<=x<=x1 and y0<=y<=y1")
parser.add_argument("--proj", nargs=1, metavar=("projname"), help="transform vertex array to projected system.\
 projname: name of the (projected) Coordinate Reference System (CRS) (e.g. EPSG:32646 for UTM46N)")
args = parser.parse_args()

if not args.objectname:
    base = os.path.basename(args.input_file)
    args.objectname = os.path.splitext(base)[0]
else:
    args.objectname = args.objectname[0]

dataxyz = np.loadtxt(args.input_file)
nvertex = np.shape(dataxyz)[0]

if args.crop:
    if args.NX:
        print("uncompatible inputs")
        exit()
    print("croping the surface")
    x0c, x1c, y0c, y1c  = [float(val) for val in args.crop]
    indexes = np.where((dataxyz[:, 0] >= x0c) & (dataxyz[:, 0] <= x1c) & (dataxyz[:, 1] >= y0c) & (dataxyz[:, 1] <= y1c))
    dataxyz = dataxyz[indexes[0], :]
    nvertex = np.shape(dataxyz)[0]

if not args.NX:
    print("NX not defined: trying to guess it...")
    rowdiff = dataxyz[1, :] - dataxyz[0, :]
    ix = -1
    ids = np.where(abs(rowdiff) < 1e-16)[0]
    if len(ids) == 1:
        # only one column starts with constant values
        ix = ids[0]
        for i in range(1, nvertex):
            if abs(dataxyz[i, ix] - dataxyz[i - 1, ix]) > 1e-16:
                NX = i
                assert nvertex % NX == 0, "nvertex%%NX!=0 nvertex/NX = %f" % (float(nvertex) / NX)
                NY = nvertex // NX
                print("NX,NY = %d,%d" % (NX, NY))
                break
    elif len(ids) > 1:
        print("2 columns starts with constant values")
        nx = []
        # find other dimension
        for ix in range(0, 3):
            if ix not in ids:
                iy = ix
        for ix in ids:
            for i in range(1, nvertex):
                if abs(dataxyz[i, ix] - dataxyz[i - 1, ix]) > 1e-16:
                    break
            if abs(dataxyz[0, iy] - dataxyz[0 + i, iy]) > 1e-16:
                nx.append(1e20)
            else:
                nx.append(i)
        NX = int(min(nx))

        if NX == 1e10:
            print("unable to guess NX and NY")
            exit()

        assert nvertex % NX == 0, "nvertex%%NX!=0 nvertex/NX = %f" % (float(nvertex) / NX)
        NY = int(nvertex / NX)
        print("NX,NY = %d,%d" % (NX, NY))
    else:
        print("unable to guess NX and NY")
        exit()
else:
    print("using user defined NX")
    NX = int(args.NX[0])
    assert nvertex % NX == 0, "nvertex%%NX!=0 nvertex/NX = %f" % (float(nvertex) / NX)
    NY = int(nvertex / NX)

if args.hole:
    print("a hole will be left in the surface")
    x0hole, x1hole, y0hole, y1hole  = [float(val) for val in args.hole]
    print("hole coordinates %f %f %f %f" % (x0hole, x1hole, y0hole, y1hole))

if args.subsample:
    onesample_every = int(args.subsample[0])
    print("subsampling : 1/%d" % onesample_every)
else:
    onesample_every = 1

dataxyz = dataxyz.reshape((NY, NX, 3))

triangles = []

dataxyz = dataxyz[::onesample_every, ::onesample_every, :]
NX = np.shape(dataxyz)[0]
NY = np.shape(dataxyz)[1]
nvertex = NX * NY


for j in range(NY - 1):
    for i in range(1, NX):
        write_triangle = True
        if args.hole:
            for ij in [[i - 1, j - 1], [i - 1, j], [i, j - 1], [i, j]]:
                if ((dataxyz[ij[0], ij[1], 0] > x0hole) & (dataxyz[ij[0], ij[1], 0] < x1hole)) & ((dataxyz[ij[0], ij[1], 1] > y0hole) & (dataxyz[ij[0], ij[1], 1] < y1hole)):
                    write_triangle = False
        if write_triangle:
            triangles.append([i + j * NX, i + 1 + j * NX, i + 1 + (j + 1) * NX])
            triangles.append([i + j * NX, i + 1 + (j + 1) * NX, i + (j + 1) * NX])

if args.proj:
    print("Projecting the nodes coordinates")
    from pyproj import Transformer
    transformer = Transformer.from_crs("epsg:4326", args.proj[0], always_xy=True)
else:
    print("no projection carried out")


fout = open(args.output_file, "w")
### WRITE THE GOCAD TS FILE
fout.write("GOCAD TSURF 1\nHEADER {\nname:" + args.objectname + "\n}\nTRIANGLES\n")
for j in range(0, NY):
    for i in range(0, NX):
        if args.proj:
            xyz = transformer.transform(dataxyz[i, j, 0], dataxyz[i, j, 1])
            fout.write(f"VRTX {i + j * NX + 1} {xyz[0]:.10e} {xyz[1]:.10e} {dataxyz[i, j, 2]:.10e}\n")
        else:
            fout.write("VRTX %d %f %f %f\n" % (i + j * NX + 1, dataxyz[i, j, 0], dataxyz[i, j, 1], dataxyz[i, j, 2]))
for tr in triangles:
    fout.write(f"TRGL {tr[0]} {tr[1]} {tr[2]}\n")
fout.write("END")

fout.close()
