#!/usr/bin/env python3
import pygmsh
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="generate a meshed box with pygmsh")
parser.add_argument("output_file", help="filename of output file")

parser.add_argument(
    "--proj",
    nargs=1,
    metavar=("projname"),
    help="transform vertex array to projected system. projname is the name of the (projected)\
Coordinate Reference System (CRS) (e.g. EPSG:32646 for UTM46N)",
)

parser.add_argument(
    "--rangeFromTopo",
    nargs=1,
    metavar=("filename"),
    help="estimate range from topography file.",
)

parser.add_argument(
    "--hdim",
    nargs=4,
    metavar=("min_lon", "max_lon", "min_lat", "max_lat"),
    help="horizontal dimension of box",
    type=float,
)
parser.add_argument(
    "--zdim",
    nargs=2,
    metavar=("zmin", "zmax"),
    help="vertical coordinates of box",
    type=float,
    required=True,
)

parser.add_argument(
    "--meshSize",
    nargs=1,
    metavar=("size"),
    help="2d mesh size",
    type=float,
    default=[10e3],
)

parser.add_argument(
    "--shrink",
    nargs=1,
    metavar=("shrink_factor"),
    help="shrink box dimension to account for the effect of projection",
    type=float,
    default=[1.0],
)

args = parser.parse_args()


def generate_projected_lon_lat_range(
    proj_string, min_longitude, max_longitude, min_latitude, max_latitude
):
    from pyproj import Transformer

    transformer = Transformer.from_crs("epsg:4326", proj_string, always_xy=True)
    central_lon = 0.5 * (min_longitude + max_longitude)
    central_lat = 0.5 * (min_latitude + max_latitude)
    __, ymin = transformer.transform(central_lon, min_latitude)
    __, ymax = transformer.transform(central_lon, max_latitude)
    xmin, __ = transformer.transform(min_longitude, central_lat)
    xmax, __ = transformer.transform(max_longitude, central_lat)
    return xmin, xmax, ymin, ymax


def read_lat_lon_range_from_netcdf(fname):
    from netCDF4 import Dataset

    with Dataset(fname, "r") as fh:
        lon = fh.variables["lon"][:]
        lat = fh.variables["lat"][:]
    return np.amin(lon), np.amax(lon), np.amin(lat), np.amax(lat)


if args.rangeFromTopo:
    if args.hdim:
        print("ignoring hdim arguments")
    min_lon, max_lon, min_lat, max_lat = read_lat_lon_range_from_netcdf(
        args.rangeFromTopo[0]
    )
elif args.hdim:
    min_lon, max_lon, min_lat, max_lat = args.hdim
else:
    raise (
        "Cannot infer box dimension: neither rangeFromTopo or hdim argument where given"
    )
if args.proj:
    xmin, xmax, ymin, ymax = generate_projected_lon_lat_range(
        args.proj[0], min_lon, max_lon, min_lat, max_lat
    )
else:
    xmin, xmax, ymin, ymax = min_lon, max_lon, min_lat, max_lat

with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_max = args.meshSize[0]
    zmin = args.zdim[0]
    xwidth = xmax - xmin
    ywidth = ymax - ymin
    height = args.zdim[1] - zmin
    box_origin = [
        xmin + 0.5 * (1 - args.shrink[0]) * xwidth,
        ymin + 0.5 * (1 - args.shrink[0]) * ywidth,
        args.zdim[0],
    ]
    box_dimension = [args.shrink[0] * xwidth, args.shrink[0] * ywidth, height]
    geom.add_box(box_origin, box_dimension, args.meshSize[0])
    mesh = geom.generate_mesh(dim=2)
    mesh.write(args.output_file)
