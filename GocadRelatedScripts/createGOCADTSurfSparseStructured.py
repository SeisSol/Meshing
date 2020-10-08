import numpy as np

# Author: Thomas Ulrich, LMU
# create surface from a structured lat/lon/depth points set
# structured in the sense: the file should have a constant step in x and y. Some values of the grid can be missing.

# parsing python arguments
import argparse
import os

parser = argparse.ArgumentParser(description="create surface from a structured grid of nodes, constant step in x and y, some nodes of the grid can be missing")
parser.add_argument("input_file", help="x y z, one node coordinate by line")
parser.add_argument("output_file", help="gocad output file")
parser.add_argument("--objectname", nargs=1, metavar=("objectname"), default=(""), help="name of the surface in gocad")
parser.add_argument("--proj", nargs=1, metavar=("projname"), default=(""), help="string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered")
args = parser.parse_args()

if args.proj != "":
    print("Projecting the nodes coordinates")
    import pyproj

    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    if args.proj[0] != "geocent":
        sProj = args.proj[0]
        myproj = pyproj.Proj(sProj)
    else:
        myproj = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")
else:
    print("no projection carried out")

if args.objectname == "":
    base = os.path.basename(args.input_file)
    args.objectname = os.path.splitext(base)[0]
else:
    args.objectname = args.objectname[0]

inputfn = args.input_file

dataxyz = np.loadtxt(inputfn)
nvertex = np.shape(dataxyz)[0]

x0 = set(dataxyz[:, 0])
x0 = sorted(x0)
dx = x0[1] - x0[0]
x = np.arange(x0[0], x0[-1] + dx, dx)

y0 = set(dataxyz[:, 1])
y0 = sorted(y0)
dy = y0[1] - y0[0]
y = np.arange(y0[0], y0[-1] + dy, dy)

nx = np.shape(x)[0]
ny = np.shape(y)[0]
xx, yy = np.meshgrid(x, y)
z = np.zeros((ny, nx))
z.fill(np.nan)

for i in range(nvertex):
    ix = int(round((dataxyz[i, 0] - x0[0]) / dx))
    iy = int(round((dataxyz[i, 1] - y0[0]) / dy))
    z[iy, ix] = dataxyz[i, 2]
triangles = []
for j in range(ny - 1):
    for i in range(nx - 1):
        if not np.any(np.isnan([z[j, i], z[j, i + 1], z[j + 1, i + 1]])):
            triangles.append([1 + i + j * nx, 1 + i + 1 + j * nx, 1 + i + 1 + (j + 1) * nx])
        if not np.any(np.isnan([z[j, i], z[j + 1, i + 1], z[j + 1, i]])):
            triangles.append([1 + i + j * nx, 1 + i + 1 + (j + 1) * nx, 1 + i + (j + 1) * nx])


fout = open(args.output_file, "w")
### WRITE THE GOCAD TS FILE
fout.write("GOCAD TSURF 1\nHEADER {\nname:" + args.objectname + "\n}\nTRIANGLES\n")
for j in range(0, ny):
    for i in range(0, nx):
        if not np.isnan(z[j, i]):
            if args.proj != "":
                xyz = pyproj.transform(lla, myproj, xx[j, i], yy[j, i], 1e3 * z[j, i], radians=False)
                fout.write("VRTX " + str(i + j * nx + 1) + " %.10e %.10e %.10e\n" % tuple(xyz))
            else:
                fout.write("VRTX %d %f %f %f\n" % (i + j * nx + 1, xx[j, i], yy[j, i], z[j, i]))
for tr in triangles:
    fout.write("TRGL %d %d %d\n" % (tr[0], tr[1], tr[2]))
fout.write("END")
fout.close()
