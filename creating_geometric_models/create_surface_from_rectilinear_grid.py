import numpy as np
from netCDF4 import Dataset
import sys
import argparse
import os
from Face import Face


class Grid:
    def __init__(self, fname, downsample):

        basename, ext = os.path.splitext(fname)
        # is there nan in the dataset (e.g. slab2 data?)
        self.is_sparse = False

        if ext == ".nc":
            self.read_netcdf(fname, downsample)
        elif ext in [".esri", ".asc"]:
            print("Warning might not work (not tested)")
            self.read_ascii_esri(fname, downsample)
        else:
            raise ValueError("format not supported", ext)
        self.compute_nx_ny()

    def read_ascii_esri(self, fname, downsample):
        "read ESRI ASCII Raster format, see e.g."
        "https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm"
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
            self.compute_nx_ny()
            print("done reading")

    def read_netcdf(self, fname, downsample):
        "read netcdf Raster file"
        fh = Dataset(fname, mode="r")
        self.keys = fh.variables.keys()
        xvar = self.determine_netcdf_variables(["lon", "x"])
        yvar = self.determine_netcdf_variables(["lat", "y"])
        zvar = self.determine_netcdf_variables(["elevation", "z", "Band1"])
        self.x = fh.variables[xvar][0::downsample]
        self.y = fh.variables[yvar][0::downsample]
        self.z = fh.variables[zvar][0::downsample, 0::downsample]
        if np.ma.is_masked(self.z):
            self.z = np.ma.filled(self.z, float("nan")) * 1e3
            self.is_sparse = True
        self.compute_nx_ny()

    def compute_nx_ny(self):
        self.nx = self.x.shape[0]
        self.ny = self.y.shape[0]
        print(self.nx, self.ny)

    def determine_netcdf_variables(self, l_possible_names):
        "search within the netcdf for variable name from a possible list of candidates"
        for name in l_possible_names:
            if name in self.keys:
                return name
        raise ("could not determine netcdf variable")

    def crop(self, argCrop):
        "crop the grid, by removing nodes out of the given latitude and longitude range"
        if argCrop:
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
            self.z = self.z.reshape(self.ny, self.nx)

    def generate_vertex(self):
        vertex = np.zeros((self.nx * self.ny, 3))
        xv, yv = np.meshgrid(self.x, self.y)
        vertex[:, 0] = xv.flatten()
        vertex[:, 1] = yv.flatten()
        vertex[:, 2] = self.z.flatten()
        self.vertex = vertex

    def proj_vertex(self, sProj):
        "project the node coordinate array"
        from pyproj import Transformer

        transformer = Transformer.from_crs("epsg:4326", args.proj[0], always_xy=True)
        print("projecting the node coordinates")
        self.vertex[:, 0], self.vertex[:, 1] = transformer.transform(self.vertex[:, 0], self.vertex[:, 1])
        print(self.vertex)
        print("done projecting")

    def translate(self):
        self.vertex[:, 0] = self.vertex[:, 0] + args.translate[0]
        self.vertex[:, 1] = self.vertex[:, 1] + args.translate[1]

    def generate_connect(self):
        if not self.is_sparse:
            self.generate_connect_full()
        else:
            self.generate_connect_sparse()

    def generate_connect_full(self):
        "triangulate the structured grid of vertex."
        "place the diagonal perpendicular to the max height gradient"
        ntriangles = 2 * (self.nx - 1) * (self.ny - 1)
        connect = np.zeros((ntriangles, 3), dtype=int)
        k = 0
        for j in range(self.ny - 1):
            for i in range(self.nx - 1):
                # edge perpendicular to the max gradient
                dz_diag1 = abs(self.vertex[i + j * self.nx, 2] - self.vertex[i + 1 + (j + 1) * self.nx, 2])
                dz_diag2 = abs(self.vertex[i + (j + 1) * self.nx, 2] - self.vertex[i + 1 + j * self.nx, 2])
                if dz_diag1 > dz_diag2:
                    connect[k, :] = [i + j * self.nx, i + 1 + j * self.nx, i + (j + 1) * self.nx]
                    connect[k + 1, :] = [i + 1 + j * self.nx, i + (j + 1) * self.nx, i + 1 + (j + 1) * self.nx]
                else:
                    connect[k, :] = [i + j * self.nx, i + 1 + j * self.nx, i + 1 + (j + 1) * self.nx]
                    connect[k + 1, :] = [i + j * self.nx, i + (j + 1) * self.nx, i + 1 + (j + 1) * self.nx]
                k = k + 2
        self.connect = connect

    def generate_connect_sparse(self):
        "triangulate the structured grid of vertex."
        "do not triangulate vertex with nan z value (used for reading Slab2.0 data)"
        triangles = []
        for j in range(self.ny - 1):
            for i in range(self.nx - 1):
                if not np.any(np.isnan(self.z[j : j + 2, i : i + 2])):
                    "place the diagonal perpendicular to the max height gradient"
                    dz_diag1 = abs(self.z[j, i] - self.z[j + 1, i + 1])
                    dz_diag2 = abs(self.z[j + 1, i] - self.z[j, i + 1])
                    if dz_diag1 > dz_diag2:
                        triangles.append([1 + i + j * self.nx, 1 + i + 1 + j * self.nx, 1 + i + (j + 1) * self.nx])
                        triangles.append([1 + i + 1 + j * self.nx, 1 + i + 1 + (j + 1) * self.nx, 1 + i + (j + 1) * self.nx])
                    else:
                        triangles.append([1 + i + j * self.nx, 1 + i + 1 + j * self.nx, 1 + i + 1 + (j + 1) * self.nx])
                        triangles.append([1 + i + j * self.nx, 1 + i + 1 + (j + 1) * self.nx, 1 + i + (j + 1) * self.nx])
                else:
                    if not np.any(np.isnan([self.z[j, i], self.z[j, i + 1], self.z[j + 1, i + 1]])):
                        triangles.append([1 + i + j * self.nx, 1 + i + 1 + j * self.nx, 1 + i + 1 + (j + 1) * self.nx])
                    elif not np.any(np.isnan([self.z[j, i], self.z[j + 1, i + 1], self.z[j + 1, i]])):
                        triangles.append([1 + i + j * self.nx, 1 + i + 1 + (j + 1) * self.nx, 1 + i + (j + 1) * self.nx])
                    elif not np.any(np.isnan([self.z[j, i], self.z[j, i + 1], self.z[j + 1, i]])):
                        triangles.append([1 + i + j * self.nx, 1 + i + 1 + j * self.nx, 1 + i + (j + 1) * self.nx])
                    elif not np.any(np.isnan([self.z[j, i + 1], self.z[j + 1, i + 1], self.z[j + 1, i]])):
                        triangles.append([1 + i + 1 + j * self.nx, 1 + i + 1 + (j + 1) * self.nx, 1 + i + (j + 1) * self.nx])
        self.connect = np.array(triangles)
        self.remove_nan_generate_vid_lookup()

    def remove_nan_generate_vid_lookup(self):
        " Generate vertex lookup array for reindexing"
        nvertex = self.vertex.shape[0]
        idv = np.linspace(0, nvertex - 1, nvertex, dtype=int)
        valid = ~np.isnan(self.vertex)[:, 2]
        # remove all entries with nan
        self.vertex = self.vertex[valid, :]
        # id of valid entries
        id_valid = idv[valid] + 1
        nvalid = id_valid.shape[0]
        # Fill in vertex lookup array for reindexing
        self.vid_lookup = {}
        for i in range(nvalid):
            self.vid_lookup[id_valid[i]] = i

    def isolate_hole(self, argHole):
        "tag a rectangular region within the grid."
        "solid_id=1 in the tagged region, 0 else"
        nconnect = self.connect.shape[0]
        self.solid_id = np.zeros(nconnect, dtype=int)
        if argHole:
            x0, x1, y0, y1 = argHole
            print("tagging hole...")
            for k in range(nconnect):
                coords = self.vertex[self.connect[k, :], 0:2]
                xmin = min(coords[:, 0])
                xmax = max(coords[:, 0])
                ymin = min(coords[:, 1])
                ymax = max(coords[:, 1])
                if ((xmin > x0) & (xmax < x1)) & ((ymin > y0) & (ymax < y1)):
                    self.solid_id[k] = 1
                else:
                    self.solid_id[k] = 0
            if max(self.solid_id) == 0:
                raise ValueError("no hole was tagged")
            print("done tagging hole")



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

structured_grid.translate()

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
