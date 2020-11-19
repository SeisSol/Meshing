from collections import defaultdict
import numpy as np


class Face:
    def __init__(self, connect):

        self.connect = connect
        self.local_vid_lookup = {}
        # for i, ivx in enumerate(self.id_vertex):
        #    self.local_vid_lookup[ivx] = i
        self.ntriangles = self.connect.shape[0]

    def compute_id_vertex(self):
        "id_vertex is an array with the id of all vertex forming the face"
        self.id_vertex = np.unique(self.connect.flatten()) + 1

    def intersect(self, Face2):
        self.compute_id_vertex()
        Face2.compute_id_vertex()
        "return the common nodes between self and Face2"
        return np.intersect1d(self.id_vertex, Face2.id_vertex)

    def reindex(self, vid_lookup):
        print("reindexing triangles...")
        for itr in range(self.ntriangles):
            for iv in range(3):
                self.connect[itr, iv] = vid_lookup[self.connect[itr, iv]]

    def writeTs(self, fname, vertex, write_full_vertex_array=True, append=False):
        """ output face as a *.ts file
             vertex: nodes coordinates array
        """
        if write_full_vertex_array:
            vertex_id_2_write = range(1, len(vertex) + 1)
        else:
            self.compute_id_vertex()
            vertex_id_2_write = self.id_vertex

        mode = "a" if append else "w"
        with open(fname, mode) as fout:
            fout.write("GOCAD TSURF 1\nHEADER {\nname:%s\nborder: true\nmesh: false\n*border*bstone: true\n}\nTFACE\n" % (fname))
            for ivx in vertex_id_2_write:
                fout.write("VRTX %s %s %s %s\n" % (ivx, vertex[ivx - 1, 0], vertex[ivx - 1, 1], vertex[ivx - 1, 2]))

            for i in range(self.ntriangles):
                fout.write("TRGL %d %d %d\n" % (self.connect[i, 0] + 1, self.connect[i, 1] + 1, self.connect[i, 2] + 1))
            fout.write("END\n")

    def computeNormal(self, vertex):
        " compute efficiently the normals "
        normal = np.cross(vertex[self.connect[:, 1], :] - vertex[self.connect[:, 0], :], vertex[self.connect[:, 2], :] - vertex[self.connect[:, 0], :])
        norm = np.apply_along_axis(np.linalg.norm, 1, normal)
        self.normal = normal / norm.reshape((self.ntriangles, 1))

    def writeStl(self, fname, vertex, append=False):
        self.computeNormal(vertex)
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

    def writebStl(self, fname, vertex):
        import struct

        fout = open(fname, "wb")
        fout.seek(80)
        fout.write(struct.pack("<L", self.ntriangles))
        self.computeNormal(vertex)
        for k in range(self.ntriangles):
            fout.write(struct.pack("<3f", *self.normal[k, :]))
            for i in range(0, 3):
                fout.write(struct.pack("<3f", *vertex[self.connect[k, i], :]))
            fout.write(struct.pack("<H", 0))

    def write(self, fname, vertex, write_full_vertex_array=True, append=False):
        import os

        basename, ext = os.path.splitext(fname)
        if ext == ".ts":
            self.writeTs(fname, vertex, write_full_vertex_array, append)
        elif ext == ".stl":
            self.writeStl(fname, vertex, append)
        elif ext == ".bstl":
            self.writebStl(fname, vertex)
        else:
            raise ValueError("format not supported", ext)
        print("done writing " + fname)
