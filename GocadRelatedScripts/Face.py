from collections import defaultdict
import numpy as np


class Face:
    def __init__(self, vertex, connect):

        self.connect = connect
        self.vertex = vertex
        self.local_vid_lookup = {}
        self.ntriangles = self.connect.shape[0]

    @classmethod
    def from_ts(self, fid):
        surface_name = "undefined"
        xyzl = []
        vid = {}
        trl = []
        prev_vert = -1
        ivertex = 0
        line = fid.readline()
        if not line:
            return None
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
            vertex = np.asarray(xyzl)
            connect = np.asarray(trl)
            myFace = Face(vertex, connect)
            myFace.reindex(vid)
            return myFace

    def setup_proj_objects(self, sProj):
        import pyproj

        lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
        if sProj != "geocent":
            myproj = pyproj.Proj(sProj)
        else:
            myproj = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")
        return lla, myproj

    def proj(self, sProj):
        "project the node coordinate array"
        import pyproj

        lla, myproj = self.setup_proj_objects(sProj)
        print("projecting the node coordinates")
        self.vertex[:, 0], self.vertex[:, 1], self.vertex[:, 2] = pyproj.transform(lla, myproj, self.vertex[:, 0], self.vertex[:, 1], self.vertex[:, 2], radians=False)
        print(self.vertex)
        print("done projecting")

    def convert_projected_to_latlon(self, sProj):
        "Convert the already projected node coordinate array to lat/lon"
        import pyproj

        lla, myproj = self.setup_proj_objects(sProj)
        print("convert the node coordinates to lat/lon")
        self.vertex[:, 0], self.vertex[:, 1], self.vertex[:, 2] = pyproj.transform(myproj, lla, self.vertex[:, 0], self.vertex[:, 1], self.vertex[:, 2], radians=False)
        print(self.vertex)
        print("done converting")

    def scale_vertex(self, scale):
        "scale vertex array"
        self.vertex = self.vertex * np.fill_diag(np.zeros((3, 3)), scale)

    def translate_vertex(self, translation):
        "translate vertex array"
        self.vertex[:, :] = self.vertex[:, :] + translation[:]

    def compute_id_vertex(self):
        "return an array with the id of all vertex forming the face"
        return np.unique(self.connect.flatten()) + 1

    def intersect(self, Face2):
        "return the common nodes between self and Face2"
        return np.intersect1d(self.compute_id_vertex(), Face2.compute_id_vertex())

    def reindex(self, vid_lookup):
        print("reindexing triangles...")
        for itr in range(self.ntriangles):
            for iv in range(3):
                self.connect[itr, iv] = vid_lookup[self.connect[itr, iv]]

    def __writeTs(self, fname, vertex, write_full_vertex_array=True, append=False):
        """ output face as a *.ts file
             vertex: nodes coordinates array
        """
        if write_full_vertex_array:
            vertex_id_2_write = range(1, len(vertex) + 1)
        else:
            vertex_id_2_write = self.compute_id_vertex()

        mode = "a" if append else "w"
        with open(fname, mode) as fout:
            fout.write("GOCAD TSURF 1\nHEADER {\nname:%s\nborder: true\nmesh: false\n*border*bstone: true\n}\nTFACE\n" % (fname))
            for ivx in vertex_id_2_write:
                fout.write("VRTX %s %s %s %s\n" % (ivx, vertex[ivx - 1, 0], vertex[ivx - 1, 1], vertex[ivx - 1, 2]))

            for i in range(self.ntriangles):
                fout.write("TRGL %d %d %d\n" % (self.connect[i, 0] + 1, self.connect[i, 1] + 1, self.connect[i, 2] + 1))
            fout.write("END\n")

    def __computeNormal(self, vertex):
        " compute efficiently the normals "
        normal = np.cross(vertex[self.connect[:, 1], :] - vertex[self.connect[:, 0], :], vertex[self.connect[:, 2], :] - vertex[self.connect[:, 0], :])
        norm = np.linalg.norm(normal, axis=1)
        self.normal = np.divide(normal, norm[:, None])

    def __writeStl(self, fname, vertex, append=False):
        self.__computeNormal(vertex)
        mode = "a" if append else "w"
        with open(fname, mode) as fout:
            fout.write(f"solid {fname}\n")
            for k in range(self.ntriangles):
                fout.write("facet normal %e %e %e\n" % tuple(self.normal[k, :]))
                fout.write("outer loop\n")
                for i in range(0, 3):
                    fout.write("vertex %.10e %.10e %.10e\n" % tuple(vertex[self.connect[k, i], :]))
                fout.write("endloop\n")
                fout.write("endfacet\n")
            fout.write(f"endsolid {fname}\n")

    def __writebStl(self, fname, vertex):
        import struct

        fout = open(fname, "wb")
        fout.seek(80)
        fout.write(struct.pack("<L", self.ntriangles))
        self.__computeNormal(vertex)
        for k in range(self.ntriangles):
            fout.write(struct.pack("<3f", *self.normal[k, :]))
            for i in range(0, 3):
                fout.write(struct.pack("<3f", *vertex[self.connect[k, i], :]))
            fout.write(struct.pack("<H", 0))

    def write(self, fname, vertex=np.empty(0), write_full_vertex_array=True, append=False):
        import os

        if vertex.shape[0] == 0:
            vertex = self.vertex

        basename, ext = os.path.splitext(fname)
        if ext == ".ts":
            self.__writeTs(fname, vertex, write_full_vertex_array, append)
        elif ext == ".stl":
            self.__writeStl(fname, vertex, append)
        elif ext == ".bstl":
            self.__writebStl(fname, vertex)
        else:
            raise ValueError("format not supported", ext)
        print("done writing " + fname)
