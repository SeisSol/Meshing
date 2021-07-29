import argparse

parser = argparse.ArgumentParser(description="project ts file")
parser.add_argument("ts_file", help="ts filename")
parser.add_argument("--proj", nargs=1, metavar=("projname"), help="transform vertex array to projected system.\
 projname: name of the (projected) Coordinate Reference System (CRS) (e.g. EPSG:32646 for UTM46N)")
args = parser.parse_args()

# set projection

from pyproj import Transformer
transformer = Transformer.from_crs("epsg:4326", args.proj[0], always_xy=True)

# read Ts file
fid = open(args.ts_file)
lines = fid.readlines()
fid.close()

foutname = args.ts_file[0:-3] + "_proj.ts"
fout = open(foutname, "w")

for line in lines:
    if line.startswith("VRTX"):
        val = [float(val) for val in line.split()[1:5]]
        xyz = transformer.transform(val[1], val[2], val[3])
        fout.write("VRTX " + str(int(val[0])) + " %.10e %.10e %.10e\n" % tuple(xyz))
    else:
        fout.write(line)
fout.close()
