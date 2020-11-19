import numpy as np
from netCDF4 import Dataset

# Author: Thomas Ulrich, LMU
# (inspired from a script by J. Klicpera)
# create surface from a structured grid of nodes

# Mesh example:

# 3 - 3
#  /
# 2 - 2
#  /
# 1 - 1

import sys

# parsing python arguments
import argparse
import os
from collections import defaultdict


class Grid:
    def __init__(self, fname, downsample):

        basename, ext = os.path.splitext(fname)
        if ext == ".nc":
            self.read_netcdf(fname, downsample)
        elif ext == ".esri":
            print("Warning might not work (not tested)")
            self.read_ascii_esri(fname, downsample)
        else:
            raise ValueError("format not supported", ext)
        self.compute_nx_ny()

    def read_ascii_esri(self, fname, downsample):
        with open(fname) as fh:
            # Read header
            NX = int(fh.readline().split()[1])
            NY = int(fh.readline().split()[1])
            xbotleft = float(fh.readline().split()[1])
            ybotleft = float(fh.readline().split()[1])
            dx = float(fh.readline().split()[1])
            fh.readline()
            # compute x and y
            self.x = np.arange(xbotleft, xbotleft + NX * dx, dx)
            self.y = np.flip(np.arange(ybotleft, ybotleft + NY * dx, dx))
            print(self.x, self.y)
            print("reading elevation")
            # using pandas rather than loadtxt because much faster
            import pandas as pd

            self.z = pd.read_csv(fh, delimiter=" ", dtype=np.float64, header=None).values[:, :-1]
            self.x = self.x[0::downsample]
            self.y = self.y[0::downsample]
            self.z = self.z[0::downsample, 0::downsample]
            print("done reading")

    def read_netcdf(self, fname, downsample):
        fh = Dataset(fname, mode="r")
        self.keys = fh.variables.keys()
        xvar = self.determine_netcdf_variables(["lon", "x"])
        yvar = self.determine_netcdf_variables(["lat", "y"])
        zvar = self.determine_netcdf_variables(["elevation", "z", "Band1"])
        self.x = fh.variables[xvar][0::downsample]
        self.y = fh.variables[yvar][0::downsample]
        self.z = fh.variables[zvar][0::downsample, 0::downsample]
        self.compute_nx_ny()

    def compute_nx_ny(self):
        self.nx = self.x.shape[0]
        self.ny = self.y.shape[0]
        print(self.nx, self.ny)

    def determine_netcdf_variables(self, l_possible_names):
        for name in l_possible_names:
            if name in self.keys:
                return name
        raise ("could not determine netcdf variable")

    def crop(self, argCrop):
        if argCrop != None:
            x0, x1, y0, y1 = argCrop
            print("crop the grid")
            print("only consider %e < lon < %e, %e < lat < %e" % (x0, x1, y0, y1))
            lon_indices = np.logical_and(self.x > x0, self.x < x1)
            lat_indices = np.logical_and(self.y > y0, self.y < y1)
            self.x = self.x[lon_indices]
            self.y = self.y[lat_indices]
            elev_indices = np.outer(lat_indices, lon_indices)
            self.z = self.z[elev_indices]
            self.compute_nx_ny()

    def generate_vertex(self):
        vertex = np.zeros((self.nx * self.ny, 3))
        xv, yv = np.meshgrid(self.x, self.y)
        vertex[:, 0] = xv.flatten()
        vertex[:, 1] = yv.flatten()
        vertex[:, 2] = self.z.flatten()
        self.vertex = vertex

    def proj_vertex(self, sProj):
        import pyproj

        lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
        if args.proj[0] != "geocent":
            myproj = pyproj.Proj(sProj)
        else:
            myproj = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")

        print("projecting the node coordinates")
        self.vertex[:, 0], self.vertex[:, 1], self.vertex[:, 2] = pyproj.transform(lla, myproj, self.vertex[:, 0], self.vertex[:, 1], self.vertex[:, 2], radians=False)
        print(self.vertex)
        print("done projecting")

    def generate_connect(self):
        ntriangles = 2 * (self.nx - 1) * (self.ny - 1)
        connect = np.zeros((ntriangles, 3), dtype=int)
        k = 0
        for j in range(self.ny - 1):
            for i in range(self.nx - 1):
                connect[k, :] = [i + j * self.nx, i + 1 + j * self.nx, i + 1 + (j + 1) * self.nx]
                connect[k + 1, :] = [i + j * self.nx, i + (j + 1) * self.nx, i + 1 + (j + 1) * self.nx]
                k = k + 2
        self.connect = connect

    def isolate_hole(self, argHole):
        nconnect = self.connect.shape[0]
        self.solid_id = np.zeros(nconnect, dtype=int)
        if argHole != None:
            x0, x1, y0, y1 = argHole
            print("tagging hole...")
            for k in range(nconnect):
                xmin, ymin = self.vertex[self.connect[k, 0], 0:2]
                xmax, ymax = self.vertex[self.connect[k, 2], 0:2]
                if ((xmin > x0) & (xmax < x1)) & ((ymin > y0) & (ymax < y1)):
                    self.solid_id[k] = 1
                    print(k)
                else:
                    self.solid_id[k] = 0
            print(max(self.solid_id))
            print("done tagging hole")


from Face import Face

parser = argparse.ArgumentParser(description="create surface from a GEBCO netcdf file")
parser.add_argument("input_file", help="GEBCO netcdf file")
parser.add_argument("output_file", help="gocad or stl output file")
parser.add_argument("--subsample", nargs=1, type=int, metavar=("onesample_every"), default=[1], help="use only one value every onesample_every in both direction")
parser.add_argument("--objectname", nargs=1, metavar=("objectname"), default=(""), help="name of the surface in gocad")
parser.add_argument("--hole", nargs=4, metavar=(("x0"), ("x1"), ("y0"), ("y1")), help="isolate a hole in surface defined by x0<=x<=x1 and y0<=y<=y1 (stl and ts output only)", type=float)
parser.add_argument("--crop", nargs=4, metavar=(("x0"), ("x1"), ("y0"), ("y1")), help="select only surfaces in x0<=x<=x1 and y0<=y<=y1", type=float)
parser.add_argument("--proj", nargs=1, metavar=("projname"), default=(""), help="string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered")
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
if args.proj != "":
    structured_grid.proj_vertex(args.proj[0])

basename, ext = os.path.splitext(args.output_file)
nsolid = max(structured_grid.solid_id) + 1

if nsolid == 1:
    myFace = Face(structured_grid.connect)
    myFace.write(f"{basename}{ext}", structured_grid.vertex)
else:
    for sid in range(nsolid):
        idtr = np.where(structured_grid.solid_id == sid)[0]
        aVid = np.unique(structured_grid.connect[idtr, :].flatten())
        myFace = Face(structured_grid.connect[idtr, :])
        myFace.write(f"{basename}{sid}{ext}", structured_grid.vertex, write_full_vertex_array=False)
