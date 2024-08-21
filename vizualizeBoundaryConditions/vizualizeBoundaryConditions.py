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
   NS=0
   SufaceId= int(args.BC)
   for faceId in range(0,4):
      boundaryFace = (boundary >> (faceId*8)) & 0xFF;
      if SufaceId==0:
         tetraId = np.where(boundaryFace>0)[0]
      elif SufaceId==-1:
         is_fault = np.invert(np.isin(boundaryFace, [0,1,5]))
         tetraId = np.where(is_fault)[0]
      else:
         tetraId = np.where(boundaryFace==SufaceId)[0]
      NS = NS + len(tetraId)
   assert(NS!=0)
   connect=np.zeros((NS,3), dtype=int)
   BC=np.zeros((1,NS))

   s_vert=np.zeros((4,3), int)
   s_vert[0,:] = [0,2,1];   s_vert[1,:] = [0,1,3];   s_vert[3,:] = [0,3,2];   s_vert[2,:] = [1,2,3];

   currentindex = 0
   for faceId in range(0,4):
      boundaryFace = (boundary >> (faceId*8)) & 0xFF;
      if SufaceId==0:
         tetraId = np.where(boundaryFace>0)[0]
      elif SufaceId==-1:
         is_fault = np.invert(np.isin(boundaryFace, [0,1,5]))
         tetraId = np.where(is_fault)[0]
      else:
         tetraId = np.where(boundaryFace==SufaceId)[0]
      for idBound in range(0, len(tetraId)):
          trinodes = tetra[ tetraId[idBound], s_vert[faceId,:]]
          connect[currentindex,:] = trinodes
          BC[0,currentindex] = boundaryFace[tetraId[idBound]]
          currentindex = currentindex +1
   print(np.unique(BC))
   aDataName = ['BC']
   prefix, ext = os.path.splitext(os.path.basename(args.filename))
   #xyz, connect, BC = remove_duplicates(xyz, connect, BC, tol=1e-4)
   sxw.write(
        f"{prefix}_bc{SufaceId}",
        xyz,
        connect,
        {"BC": BC},
        {},
        reduce_precision=True,
        backend="hdf5",
    )

parser = argparse.ArgumentParser(
    description="Read a PUML mesh and create a xdmf/h5 file containing the surface boundary mesh"
)
parser.add_argument(
    "filename",
    help="PUML mesh",
)
parser.add_argument("BC", help="1: free surface, 3:dynamic rupture 5:absorbing 0:all")
args = parser.parse_args()

ReadHdf5PosixForBoundaryPlotting(args.filename)
