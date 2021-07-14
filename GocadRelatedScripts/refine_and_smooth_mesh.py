from Face import Face
import argparse
import os
import trimesh

parser = argparse.ArgumentParser(description="refine and smooth")
parser.add_argument("input_file", help="surface in ts file format")
parser.add_argument("--N",
                    nargs=1,
                    metavar=("N"),
                    type=int,
                    default=0,
                    help="number of initial refine steps")
parser.add_argument("--P",
                    nargs=1,
                    metavar=("P"),
                    type=int,
                    default=1,
                    help="number of and refine/smoothing steps")
args = parser.parse_args()

fid = open(args.input_file)
myFace = Face.from_ts(fid)

a = trimesh.Trimesh(vertices=myFace.vertex, faces=myFace.connect)

for i in range(args.N[0]):
    a = a.subdivide()
for i in range(args.P[0]):
    a = a.subdivide()
    a = trimesh.smoothing.filter_taubin(a)

myFace = Face(a.vertices, a.faces)
basename, ext = os.path.splitext(args.input_file)
myFace.write(f"{basename}_refined_smooth_{args.N[0]}{ext}")
