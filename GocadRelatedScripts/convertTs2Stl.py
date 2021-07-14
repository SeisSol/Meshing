import numpy as np
import argparse
import sys
import os
from Face import Face

parser = argparse.ArgumentParser(description="convert ts file to stl")
parser.add_argument("ts_file", help="ts filename")
parser.add_argument(
    "--proj",
    nargs=1,
    metavar=("projname"),
    default=None,
    help=
    "string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered"
)
parser.add_argument("--tokm",
                    dest="tokm",
                    action="store_true",
                    help="convert coordinates to km")
parser.add_argument("--translate",
                    nargs=2,
                    metavar=("x0", "y0"),
                    default=([0, 0]),
                    help="translates all nodes by (x0,y0)",
                    type=float)
args = parser.parse_args()

idface = 0
basename, ext = os.path.splitext(args.ts_file)

with open(args.ts_file) as fid:
    while True:
        myFace = Face.from_ts(fid)
        if myFace == None:
            break
        if args.proj != None:
            myFace.proj(args.proj[0])
        if args.tokm:
            myFace.scale_vertex([0.001, 0.001, 0.001])
        if args.translate:
            myFace.translate_vertex([args.translate[0], args.translate[1], 0])
        myFace.write(f"{basename}.stl", append=(idface != 0))
        idface += 1
