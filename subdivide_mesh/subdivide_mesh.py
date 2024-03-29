#!/usr/bin/env python3
import seissolxdmf as sx
import seissolxdmfwriter as sw
import numpy as np
import re
import argparse
from numba import jit, njit, prange
import os

parse_version = re.compile("(\d+)\.(\d+)\.(\d+)")
version = [int(v) for v in parse_version.match(sw.__version__).groups()]
if (version[0] == 0) and (version[1] < 2):
    print(
        f"Need at least seissolxdmfwriter version 0.2.0, but found {version[0]}.{version[1]}.{version[2]}"
    )
    quit()


def read_xdmf_mesh(filename):
    mesh = sx.seissolxdmf(filename)
    geom = mesh.ReadGeometry()
    connect = mesh.ReadConnect()
    boundary = mesh.ReadData("boundary")
    group = mesh.ReadData("group")
    return geom, connect, boundary, group


def write_xdmf_mesh(filename, geom, connect, boundary, group):
    sw.write_seissol_output(
        filename, geom, connect, ["boundary", "group"], [boundary, group], 0.0, [0]
    )


@njit(cache=True)
def unpack_boundary(boundary):
    return np.array([0xFF & (boundary >> (8 * i)) for i in range(4)]).flatten()


@njit(cache=True)
def pack_boundary(boundary):
    return boundary[0] + (boundary[1] << 8) + (boundary[2] << 16) + (boundary[3] << 24)


@njit(cache=True)
def subdivide_element(
    index,
    geom,
    connect,
    boundary,
    group,
    new_geom,
    new_connect,
    new_boundary,
    new_group,
):
    i = index
    element = connect[index, :]
    vertices = geom[element]
    a = vertices[0, :]
    b = vertices[1, :]
    c = vertices[2, :]
    d = vertices[3, :]
    ab = 0.5 * (a + b)
    ac = 0.5 * (a + c)
    ad = 0.5 * (a + d)
    bc = 0.5 * (b + c)
    bd = 0.5 * (b + d)
    cd = 0.5 * (c + d)
    # new_geom[10 * i : 10 * (i + 1), :] = np.array([a, b, c, d, ab, ac, ad, bc, bd, cd])
    new_geom[10 * i + 0, :] = a
    new_geom[10 * i + 1, :] = b
    new_geom[10 * i + 2, :] = c
    new_geom[10 * i + 3, :] = d
    new_geom[10 * i + 4, :] = ab
    new_geom[10 * i + 5, :] = ac
    new_geom[10 * i + 6, :] = ad
    new_geom[10 * i + 7, :] = bc
    new_geom[10 * i + 8, :] = bd
    new_geom[10 * i + 9, :] = cd
    new_connect[8 * i : 8 * (i + 1), :] = (
        np.array(
            [
                [0, 4, 5, 6],
                [1, 4, 8, 7],
                [2, 5, 7, 9],
                [3, 6, 9, 8],
                [4, 5, 6, 8],
                [4, 5, 8, 7],
                [5, 6, 8, 9],
                [5, 7, 9, 8],
            ]
        )
        + 10 * i
    )

    original_boundary = unpack_boundary(boundary[index])
    modified_boundary = np.array(
        [
            [original_boundary[0], original_boundary[1], 0, original_boundary[3]],
            [original_boundary[1], original_boundary[0], 0, original_boundary[2]],
            [original_boundary[0], original_boundary[3], 0, original_boundary[2]],
            [original_boundary[3], original_boundary[1], 0, original_boundary[2]],
            [0, 0, 0, original_boundary[1]],
            [0, original_boundary[0], 0, 0],
            [0, original_boundary[3], 0, 0],
            [0, 0, original_boundary[2], 0],
        ]
    )

    new_boundary[8 * i : 8 * (i + 1)] = np.array(
        [pack_boundary(b) for b in modified_boundary]
    )
    new_group[8 * i : 8 * (i + 1)] = np.array([group[index]] * 8)


@njit(parallel=False, cache=True)
def do_subdivide(geom, connect, boundary, group):
    number_of_elements = connect.shape[0]
    new_geom = np.zeros((10 * number_of_elements, 3), dtype=np.float64)
    new_connect = np.zeros((8 * number_of_elements, 4), dtype=np.uint64)
    new_boundary = np.zeros((8 * number_of_elements,), dtype=np.int32)
    new_group = np.zeros((8 * number_of_elements,), dtype=np.int32)

    for i in prange(number_of_elements):
        subdivide_element(
            i,
            geom,
            connect,
            boundary,
            group,
            new_geom,
            new_connect,
            new_boundary,
            new_group,
        )

    return new_geom, new_connect, new_boundary, new_group


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read xdmf mesh and subdivide each element into eight subelements"
    )
    parser.add_argument("filename", help="xdmf mesh file")
    args = parser.parse_args()

    parse_filename = re.compile("(.+)\.xdmf")

    filename_prefix = parse_filename.match(args.filename)
    if filename_prefix == None:
        print("Need to specify an xdmf file")
        quit()
    prefix_out = os.path.basename(f"{filename_prefix.group(1)}-subdivided")

    geom, connect, boundary, group = read_xdmf_mesh(args.filename)

    new_geom, new_connect, new_boundary, new_group = do_subdivide(
        geom, connect, boundary, group
    )
    # Needs to be outside jit
    new_geom, inverse = np.unique(new_geom, return_inverse=True, axis=0)
    new_connect = inverse[new_connect]

    write_xdmf_mesh(prefix_out, new_geom, new_connect, new_boundary, new_group)
