#!/usr/bin/env python3
import seissolxdmf
import seissolxdmfwriter as sxw
import numpy as np
import argparse
import os

def remove_duplicates(geom, connect, Bc, tol=1e-4):
    """sort geom array and reindex connect array to match the new geom array"""
    import pymesh

    nv = geom.shape[0]
    geom, connect, inf = pymesh.remove_duplicated_vertices_raw(geom, connect)
    geom, connect, inf = pymesh.remove_duplicated_faces_raw(geom, connect)
    ori_face_index = inf['ori_face_index']
    BC = BC[ori_face_index]
    return geom, connect, BC


def ReadHdf5PosixForBoundaryPlotting(filename):
   sx = seissolxdmf.seissolxdmf(filename)
   xyz = sx.ReadGeometry()
   tetra = sx.ReadConnect()
   nElements=np.shape(tetra)[0]
   boundary = sx.ReadDataChunk('boundary', firstElement=0, nchunk=nElements, idt=0)

   def get_tetras_with_passing_boundary_condition(boundary,faceId, BC):
          boundaryFace = (boundary >> (faceId*8)) & 0xFF;
          if BC == "all":
             return np.where(boundaryFace>0)[0]
          elif BC == "faults":
             is_fault = np.invert(np.isin(boundaryFace, [0,1,5]))
             return np.where(is_fault)[0]
          else:
             return np.where(boundaryFace==BC)[0]

   def compute_output_mesh_size(boundary):
       NS=0
       for faceId in range(0,4):
          tetraId = get_tetras_with_passing_boundary_condition(boundary,faceId, args.BC)
          NS = NS + len(tetraId)
       assert(NS!=0)
       return NS

   NS = compute_output_mesh_size(boundary)
   connect=np.zeros((NS,3), dtype=int)
   BC=np.zeros((1,NS))

   s_vert=np.zeros((4,3), int)
   s_vert[0,:] = [0,2,1];   s_vert[1,:] = [0,1,3];   s_vert[3,:] = [0,3,2];   s_vert[2,:] = [1,2,3];

   currentindex = 0
   for faceId in range(0,4):
      boundaryFace = (boundary >> (faceId*8)) & 0xFF;
      tetraId = get_tetras_with_passing_boundary_condition(boundary,faceId, args.BC)
      for idBound in range(0, len(tetraId)):
          trinodes = tetra[ tetraId[idBound], s_vert[faceId,:]]
          connect[currentindex,:] = trinodes
          BC[0,currentindex] = boundaryFace[tetraId[idBound]]
          currentindex = currentindex +1
   print(np.unique(BC))
   prefix, ext = os.path.splitext(os.path.basename(args.filename))
   #xyz, connect, BC = remove_duplicates(xyz, connect, BC, tol=1e-4)
   sBC = "BC" if BC != "faults"  else "fault-tag"
   sxw.write(
        f"{prefix}_bc_{args.BC}",
        xyz,
        connect,
        {"BC": BC},
        {},
        reduce_precision=True,
        backend="hdf5",
    )

def custom_choice(value):
    # Allow "all" and "fault" directly
    allowed_strings = ["all", "free-surface", "absorbing", "faults"]
    if value == "free-surface":
        return 1
    elif value == "absorbing":
        return 5
    elif value in allowed_strings:
        return value
    # Try to convert the value to an integer
    try:
        ivalue = int(value)
        return ivalue
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid choice: '{value}' (choose from {allowed_strings}, or an integer)")


parser = argparse.ArgumentParser(
    description="Read a PUML mesh and create a xdmf/h5 file containing the surface boundary mesh"
)
parser.add_argument(
    "filename",
    help="PUML mesh",
)
parser.add_argument("BC", 
    type=custom_choice,
    help="boundary condition can be 'all', 'free-surface', 'faults', 'absorbing' or an integer.")

args = parser.parse_args()

ReadHdf5PosixForBoundaryPlotting(args.filename)
