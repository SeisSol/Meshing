import numpy as np
import argparse
import sys
import os
from Face import Face


parser = argparse.ArgumentParser(description="convert ts file to another supported format (e.g. stl, bstl, ts)")
parser.add_argument("ts_file", help="ts filename")
parser.add_argument("output_filename", nargs="?", help="output filname (if not used = basename.stl)")
parser.add_argument("--proj", nargs=1, metavar=("projname"), help="transform vertex array to projected system.\
projname: name of the (projected) Coordinate Reference System (CRS) (e.g. EPSG:32646 for UTM46N)")
parser.add_argument("--tokm", dest="tokm", action="store_true", help="convert coordinates to km")
parser.add_argument("--translate", nargs=2, metavar=("x0", "y0"), default=([0, 0]), help="translates all nodes by (x0,y0)", type=float)
args = parser.parse_args()


idface = 0
if not args.output_filename:
    args.output_filename = args.ts_file[0:-3] + ".stl"
basename, ext = os.path.splitext(os.path.basename(args.output_filename))

with open(args.ts_file) as fid:
    while True:
        myFace = Face.from_ts(fid)
        if not myFace:
            break
        if args.proj:
            myFace.proj(args.proj[0])
        if args.tokm:
            myFace.scale_vertex([0.001, 0.001, 0.001])
        if args.translate:
            myFace.translate_vertex([args.translate[0], args.translate[1], 0])
        if ext in [".stl", ".ts"]:
            myFace.write(f"{basename}{ext}", append=(idface != 0))
        else:
            myFace.write(f"{basename}{idface}{ext}")
        idface += 1
