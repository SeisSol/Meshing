#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from Grid import Grid
from math import floor

parser = argparse.ArgumentParser(
    description="read trench file, get depth from gebco file, shift"
)
parser.add_argument("output_file", help="filename of output file")
parser.add_argument(
    "--usgs_trench_file",
    nargs=1,
    metavar=("filename"),
    help="pass to usgs trench filename. \
Available at https://github.com/usgs/slab2/blob/master/misc/trenches/trenches_usgs_2017.csv",
    required=True,
)
parser.add_argument(
    "--filter_trench",
    nargs=5,
    metavar=("min_lon", "max_lon", "min_lat", "max_lat", "trench_name"),
    help="filter trench by location and name",
)

parser.add_argument(
    "--bathymetry",
    nargs=1,
    metavar=("filename"),
    help="path to (gebco) bathymetry file.",
    required=True,
)
parser.add_argument(
    "--shift",
    nargs=2,
    metavar=("x", "y"),
    help="shift trench location new = old + (x,y) degrees",
    type=float,
    required=True,
)
parser.add_argument(
    "--plot",
    dest="plot",
    action="store_true",
    default=False,
    help="plot elevation at the trench vs latitude and longitude vs latitude",
)
args = parser.parse_args()

bgrid = Grid(args.bathymetry[0], 1)


def extract_trench(args):
    trench = pd.read_csv(args.usgs_trench_file[0])
    xmin, xmax, ymin, ymax = [float(val) for val in args.filter_trench[0:4]]
    trench_name = args.filter_trench[4]
    subtrench = trench[
        (trench["lon"].between(xmin, xmax))
        & (trench["lat"].between(ymin, ymax))
        & (trench["slab"] == trench_name)
    ]
    subtrench.loc[:, "z"] = 0.0
    return subtrench


subtrench = extract_trench(args)

# Update trench elevation with elevation from bathymetry file
lon0 = bgrid.x[0]
dlon = bgrid.x[1] - bgrid.x[0]
lat0 = bgrid.y[0]
dlat = bgrid.y[1] - bgrid.y[0]
for index, row in subtrench.iterrows():
    i = floor((row["lon"] - lon0) / dlon)
    j = floor((row["lat"] - lat0) / dlat)
    subtrench.loc[index, "z"] = bgrid.z[j, i]

if args.plot:
    # plot seafloor depth
    import matplotlib.pylab as plt

    f, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(subtrench["lat"].values, subtrench["z"].values)
    ax1.set_xlabel("latitude")
    ax1.set_ylabel("elevation along trench (m)")
    ax2.plot(subtrench["lon"].values, subtrench["lat"].values)
    ax2.set_xlabel("longitude")
    ax2.set_ylabel("latitude")
    ax2.set_aspect("equal")
    plt.show()

subtrench.loc[:]["lon"] += args.shift[0]
subtrench.loc[:]["lat"] += args.shift[1]

df = subtrench[["lon", "lat", "z"]]
np.savetxt(args.output_file, df.values, fmt="%f")
