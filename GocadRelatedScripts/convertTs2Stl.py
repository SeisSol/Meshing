import numpy as np
import argparse
import sys
import os
from Face import Face


def ReadBlockGocadTs(fid):
    surface_name = "undefined"
    xyzl = []
    vid = {}
    trl = []
    prev_vert = -1
    ivertex = 0
    line = fid.readline()
    if not line:
        return False
    title = line.split()
    # Skip everything not Surface
    if title[1].lower() != "tsurf":
        print("skipping %s" % (title[1]))
        while True:
            line = fid.readline()
            if not line:
                break
            if line.startswith("END_"):
                continue
            elif line.startswith("END"):
                break
    else:
        while True:
            line = fid.readline()
            if not line:
                break
            if line.startswith("VRTX"):
                val = [float(val) for val in line.split()[1:5]]
                newVid = int(val[0])
                vid[newVid] = ivertex
                ivertex = ivertex + 1
                xyzl.append(val[1:4])
            elif line.startswith("TRGL"):
                val = [int(val) for val in line.split()[1:4]]
                trl.append(val)
            elif line.startswith("END_"):
                continue
            elif line.startswith("name"):
                surface_name = (line.split(":")[1]).strip()
                print("now processing %s" % surface_name)
            elif line.startswith("ATOM"):
                val = [int(val) for val in line.split()[1:3]]
                vid0 = vid[val[1]]
                xyzl.append(xyzl[vid0])
                vid[val[0]] = ivertex
                ivertex = ivertex + 1
            elif line.startswith("END"):
                break
        nodes = np.asarray(xyzl)
        myFace = Face(np.asarray(trl))
        myFace.reindex(vid)
        return [nodes, surface_name, myFace]


parser = argparse.ArgumentParser(description="convert ts file to stl")
parser.add_argument("ts_file", help="ts filename")
parser.add_argument("--proj", nargs=1, metavar=("projname"), default=(""), help="string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered")
parser.add_argument("--tokm", dest="tokm", action="store_true", help="convert coordinates to km")
parser.add_argument("--translate", nargs=2, metavar=("x0", "y0"), default=([0, 0]), help="translates all nodes by (x0,y0)", type=float)
args = parser.parse_args()


# set projection
if args.proj != "":
    import mpl_toolkits.basemap.pyproj as pyproj

    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    if args.proj[0] != "geocent":
        sProj = args.proj[0]
        myproj = pyproj.Proj(sProj)
    else:
        myproj = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")

idface = 0
basename, ext = os.path.splitext(args.ts_file)

with open(args.ts_file) as fid:
    while True:
        outputs = ReadBlockGocadTs(fid)
        if outputs == False:
            break
        else:
            nodes, surface_name, myFace = outputs
        print("done reading %s, found %d nodes and %d triangles" % (surface_name, nodes.shape[0], myFace.ntriangles))
        if args.proj != "":
            print("projecting the nodes coordinates")
            xyzb = pyproj.transform(lla, myproj, nodes[:, 0], nodes[:, 1], 1e3 * nodes[:, 2], radians=False)
            nodes[:, 0] = xyzb[0]
            nodes[:, 1] = xyzb[1]
            nodes[:, 2] = xyzb[2]
        nodes[:, 0] = nodes[:, 0] + args.translate[0]
        nodes[:, 1] = nodes[:, 1] + args.translate[1]
        if args.tokm:
            nodes[:, :] = nodes[:, :] / 1e3
        myFace.write(f"{basename}.stl", nodes, append=(idface != 0))
        idface += 1
