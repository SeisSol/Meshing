import numpy as np
from math import tan, pi, pow, cos, sin
from scipy.interpolate import interp1d
import os
from datetime import datetime

import argparse

parser = argparse.ArgumentParser(description="create curved fault geometry from pl file")
parser.add_argument("filename", help="fault trace (*.pl) or ascii file (2 or 3 columns)")
parser.add_argument("dipType", type=int, help="0: constant dip, 1: depth dependant dip, described by an ascii file, 2: dip variying along the length of the trace")
parser.add_argument("dipDesc", help="dipType=0: dip value dipType=1 name of ascii file with 2 columns (depth, dip). dipType=2: idem with (relative length[0-1], dip)")
parser.add_argument("--constExtrudeDir", nargs=1, metavar=("strike"), help="dip horizontal direction defined by given strike (else changing along trace). if strike=360, extreme points of curve define the strike", type=float)
# todo: average strike instead of 2 extreme points?
parser.add_argument("--translate", nargs=2, metavar=("x0", "y0"), default=([0, 0]), help="translates all nodes by (x0,y0)", type=float)
parser.add_argument("--dd", nargs=1, metavar=("dd"), default=([1e3]), help="sampling along depth (m)", type=float)
parser.add_argument("--maxdepth", nargs=1, metavar=("maxdepth"), default=([20e3]), help="max depth (positive) of fault (m)", type=float)
parser.add_argument("--extend", nargs=1, metavar=("extend"), default=([00e3]), help="extend toward z= extend (positive)", type=float)
parser.add_argument("--first_node_ext", nargs=1, default=([0.0]), help="extend along strike trace before the first node (arg: length in km)", type=float)
parser.add_argument("--last_node_ext", nargs=1, default=([0.0]), help="extend along strike trace after the last node (arg: length in km)", type=float)
parser.add_argument("--proj", nargs=1, metavar=("projname"), default=(""), help="string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered")
parser.add_argument("--smoothingParameter", nargs=1, metavar=("smoothingParameter"), default=([1e5]), help="smoothing parameter (the bigger the smoother)", type=float)
parser.add_argument("--plotFaultTrace", dest="plotFaultTrace", action="store_true", default=False, help="plot resampled fault trace in matplotlib")
args = parser.parse_args()


def Project(nodes):
    print("Projecting the nodes coordinates")
    import pyproj

    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    if args.proj[0] != "geocent":
        sProj = args.proj[0]
        myproj = pyproj.Proj(sProj)
    else:
        myproj = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")
    xyz = pyproj.transform(lla, myproj, nodes[:, 0], nodes[:, 1], 1e3 * nodes[:, 2], radians=False)
    nodes[:, 0] = xyz[0]
    nodes[:, 1] = xyz[1]
    nodes[:, 2] = xyz[2]
    return nodes


dx = args.dd[0]


# Reading dip value
if args.dipType == 0:
    dip = float(args.dipDesc) * pi / 180.0
elif args.dipType == 1:
    print("depth dependant dip described (depth vs dip) by the file %s" % (args.dipDesc))
    depthdip = np.loadtxt(args.dipDesc)
    deptha = depthdip[:, 0]
    dipa = depthdip[:, 1] * pi / 180.0
    dipangle = interp1d(deptha, dipa, kind="linear")
elif args.dipType == 2:
    print("along strike varying dip described (relative length along strike (0-1) vs dip) by file %s" % (args.dipDesc))
    curviligndip = np.loadtxt(args.dipDesc)
    relD = curviligndip[:, 0]
    dipa = curviligndip[:, 1] * pi / 180.0
    dipangle = interp1d(relD, dipa, kind="linear")
else:
    print("dipType not in 0-2")
    print(args.dipType)
    exit()


# Reading fault trace
bn = os.path.basename(args.filename)
ext = bn.split(".")[1]
# Trace described by a pl file
if ext == "pl":
    nodes = []
    fid = open(args.filename)
    lines = fid.readlines()
    fid.close()
    for li in lines:
        if li.startswith("VRTX"):
            lli = li.split()
            nodes.append([float(lli[2]), float(lli[3]), float(lli[4])])
    nodes = np.asarray(nodes)
# Trace described by an Ascii file
else:
    nodes = np.loadtxt(args.filename)
    ndim = np.shape(nodes)[1]
    if ndim == 2:
        nx = np.shape(nodes)[0]
        b = np.zeros((nx, 1))
        nodes = np.append(nodes, b, axis=1)

if args.proj != "":
    nodes = Project(nodes)

if args.first_node_ext[0] > 0:
    u0 = nodes[0, :] - nodes[1, :]
    u0 = u0 / np.linalg.norm(u0)
    nodes = np.vstack([nodes[0, :] + u0 * args.first_node_ext[0], nodes])

if args.last_node_ext[0] > 0:
    u0 = nodes[-1, :] - nodes[-2, :]
    u0 = u0 / np.linalg.norm(u0)
    nodes = np.vstack([nodes, nodes[-1, :] + u0 * args.last_node_ext[0]])


nodes[:, 0] = nodes[:, 0] + args.translate[0]
nodes[:, 1] = nodes[:, 1] + args.translate[1]

# Compute Fault length and nx
diff = nodes[1:, :] - nodes[0:-1, :]
faultlength = np.sum(np.sqrt(np.square(diff[:, 0]) + np.square(diff[:, 1])))
print("faultlength = %.2f km" % (faultlength / 1e3))
nx = int(faultlength / dx)

# smooth and resample fault trace
from scipy.interpolate import splprep, splev

tck, u = splprep([nodes[:, 0], nodes[:, 1], nodes[:, 2]], s=args.smoothingParameter[0])
unew = np.linspace(0, 1, nx)
new_points = splev(unew, tck)
nNewNodes = np.shape(new_points[0])[0]
nodes = np.zeros((nNewNodes, 3))
nodes[:, 0] = new_points[0]
nodes[:, 1] = new_points[1]
nodes[:, 2] = new_points[2]

# Plot smooth fault trace over picked trace
if args.plotFaultTrace:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(nodes[:, 0], nodes[:, 1], "ro")
    plt.axis("equal")
    ax.plot(new_points[0], new_points[1], "rx-")
    plt.show()


# Compute x,y abscisse coordinates of each nodes
diff = nodes[1:, :] - nodes[0:-1, :]
distBetweenNodes = np.sqrt(np.square(diff[:, 0]) + np.square(diff[:, 1]))
print("distance between nodes: 15, 50 and 85 percentiles", np.percentile(distBetweenNodes, 15), np.percentile(distBetweenNodes, 50), np.percentile(distBetweenNodes, 85))

xi = np.zeros((nNewNodes))
for i in range(1, nNewNodes):
    xi[i] = xi[i - 1] + distBetweenNodes[i - 1]


depth = -np.flip(np.arange(0, args.maxdepth[0], dx), axis=0)

nx = np.shape(nodes)[0]
nd = np.shape(depth)[0]
vertices = np.zeros((nx, nd, 3))
uz = np.array([0, 0, 1])
vertices[:, 0, :] = nodes

if args.dipType == 2:
    # we need to compute the curvilign coordinate along strike
    dist = np.linalg.norm(nodes[1:nx] - nodes[0 : nx - 1], axis=1)
    distall = np.sum(dist)
    reldist_seg = dist / distall
    reldist = np.zeros(nx)
    for i in range(1, nx):
        reldist[i] = reldist[i - 1] + reldist_seg[i - 1]
    reldist[reldist > 1] = 1.0

# v0: unit strike vector
# v1: normal to v0, in the horizontal plane
av0 = np.zeros((nx, 3))
av1 = np.zeros((nx, 3))
for i in range(0, nx):
    if args.constExtrudeDir == None:
        if i + 1 != nx:
            v0 = vertices[i + 1, 0, :] - vertices[i, 0, :]
        else:
            v0 = vertices[i, 0, :] - vertices[i - 1, 0, :]
    else:
        if args.constExtrudeDir[0] == 360:
            v0 = vertices[nx - 1, 0, :] - vertices[0, 0, :]
        else:
            strike = args.constExtrudeDir[0] * pi / 180.0
            v0 = [-cos(strike), -sin(strike), 0.0]
    v0[2] = 0
    v0 = v0 / np.linalg.norm(v0)
    av0[i, :] = v0
    av1[i, :] = np.array([-v0[1], v0[0], 0])

# Create new vertex below 0
for i in range(0, nx):
    for j in range(1, nd):
        if args.dipType == 0:
            mydip = dip
        elif args.dipType == 1:
            mydip = dipangle(depth[nd - j] + vertices[i, 0, 2])
        else:
            mydip = dipangle(reldist[i])
        ud = -(1.0 / tan(mydip) * av1[i, :] - uz)
        vertices[i, j, :] = vertices[i, j - 1, :] - dx * ud


# Create new vertex above 0
if args.extend[0] > 0:
    # Extension toward z plus
    depth = np.flip(np.arange(0, args.extend[0], dx), axis=0)
    # depth = np.arange(0,args.extend[0],dx)
    nd2 = np.shape(depth)[0]
    vertices2 = np.zeros((nx, nd2, 3))
    vertices2[:, nd2 - 1, :] = nodes
    for i in range(0, nx):
        for j in range(1, nd2):
            if args.dipType == 0:
                mydip = dip
            elif args.dipType == 1:
                mydip = dipangle(depth[nd2 - j] + vertices[i, 0, 2])
            else:
                mydip = dipangle(reldist[i])
            ud = -(1.0 / tan(mydip) * av1[i, :] - uz)
            normal = np.cross(v0, ud)
            normal = normal / np.linalg.norm(v0)
            vertices2[i, nd2 - 1 - j, :] = vertices2[i, nd2 - j, :] + dx * ud
    vertices = np.concatenate((vertices2, vertices[:, 1:, :]), axis=1)

if args.extend[0] > 0:
    # Extension toward z plus
    depth = np.flip(np.arange(0, args.extend[0], dx), axis=0)
    depth = np.concatenate((depth, -np.arange(dx, args.maxdepth[0], dx)))
else:
    depth = -np.arange(dx, args.maxdepth[0], dx)
nd = np.shape(depth)[0]

### WRITE THE GOCAD TS FILE
bn = os.path.basename(args.filename)
prefix = bn.split(".")[0]
NX = nx
NY = nd

fout = open(prefix + "0.ts", "w")
fout.write("GOCAD TSURF 1\nHEADER {\nname:" + prefix + "\n}\nTRIANGLES\n")
for j in range(0, NY):
    for i in range(0, NX):
        fout.write("VRTX " + str(i + j * NX + 1) + " %.10e %.10e %.10e\n" % tuple(vertices[i, j, :]))
for j in range(NY - 1):
    for i in range(1, NX):
        fout.write("TRGL %d %d %d\n" % (i + j * NX, i + 1 + j * NX, i + 1 + (j + 1) * NX))
        fout.write("TRGL %d %d %d\n" % (i + j * NX, i + 1 + (j + 1) * NX, i + (j + 1) * NX))
fout.write("END")

fout.close()
