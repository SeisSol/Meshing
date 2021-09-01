#!/usr/bin/env python3
import seissolxdmf
import numpy as np
import argparse
import os

def CreateXdmf(fn, nNodes, ntriangles, aDataName):
   xdmf="""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
   <Grid Name="step_000000000000" GridType="Uniform"><!-- mesh id: 0, mesh step: 0 -->
    <Topology TopologyType="Triangle" NumberOfElements="%d">
     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="%d 3">%s.h5:/mesh0/connect</DataItem>
    </Topology>
    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="%d">
     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="%d 3">%s.h5:/mesh0/geometry</DataItem>
    </Geometry>
    <Time Value="0"/>""" %(ntriangles, ntriangles, fn, nNodes, nNodes, fn)

   for dataName in aDataName:
      xdmf=xdmf + """
    <Attribute Name="%s" Center="Cell">
     <DataItem ItemType="HyperSlab" Dimensions="%d">
      <DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">0 0 1 1 1 %d</DataItem>
      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 %d">%s.h5:/mesh0/%s</DataItem>
     </DataItem>
    </Attribute>""" % (dataName, ntriangles, ntriangles, ntriangles,fn,dataName)

   xdmf=xdmf + """
   </Grid>
  </Grid>
 </Domain>
</Xdmf>
"""
   fid=open(fn+'.xdmf','w')
   fid.write(xdmf)
   fid.close()

def write_xdmfh5(fname, aDataName, xyz, connect, BC):
   nNodes=xyz.shape[0]
   ntriangles = connect.shape[0]
   CreateXdmf(fname, nNodes, ntriangles, aDataName)
   #Write h5 file
   import h5py
   h5f = h5py.File(fname+'.h5','w')
   h5f.create_dataset('mesh0/connect', (ntriangles,3), dtype='i8')
   h5f['mesh0/connect'][:,:] = connect[:,:]
   h5f.create_dataset('mesh0/geometry', xyz.shape, dtype='d')
   h5f['mesh0/geometry'][:,:] = xyz[:,:]
   for dataName in aDataName:
     hdname = "mesh0/"+dataName
     h5f.create_dataset(hdname, (1, ntriangles), dtype='d')
     h5f[hdname][0,:] = eval(dataName)[:]
   h5f.close()
   print ("done writing %s.xdmf" %fname)



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
      else:
         tetraId = np.where(boundaryFace==SufaceId)[0]
      for idBound in range(0, len(tetraId)):
          trinodes = tetra[ tetraId[idBound], s_vert[faceId,:]]
          connect[currentindex,:] = trinodes
          BC[0,currentindex] = boundaryFace[tetraId[idBound]]
          currentindex = currentindex +1

   aDataName = ['BC']
   prefix, ext = os.path.splitext(args.filename)
   fn = prefix+'_bc%d' %(SufaceId)
   write_xdmfh5(fn, aDataName, xyz, connect, BC)

parser = argparse.ArgumentParser(description='Read hdf5 mesh and create a xdmf/h5 file containing the BC surfaces')
parser.add_argument('filename', help='fault output filename (xdmf), or SeisSol netcdf (nc) or ts (Gocad)')
parser.add_argument('BC', help='1: free surface, 3:dynamic rupture 5:absorbing 0:all')
args = parser.parse_args()

ReadHdf5PosixForBoundaryPlotting(args.filename)
