#!/usr/bin/env python3
import seissolxdmf
import numpy as np
import argparse
import os

def CreateXdmf(fn, Nnodesb, nElementsb):
   xdmf="""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="puml mesh" GridType="Uniform">
   <Geometry name="geo" GeometryType="XYZ" NumberOfElements="%d">
    <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="%d 3">%s:/geometry</DataItem>
   </Geometry>
   <Topology TopologyType="Tetrahedron" NumberOfElements="%d">
    <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="%d 4">%s:/connect</DataItem>
   </Topology>
   <Attribute Name="group" Center="Cell">
    <DataItem  NumberType="Int" Precision="4" Format="HDF" Dimensions="%d">%s:/group</DataItem>
   </Attribute>
   <Attribute Name="boundary" Center="Cell">
    <DataItem NumberType="Int" Precision="4" Format="HDF" Dimensions="%d">%s:/boundary</DataItem>
   </Attribute>
  </Grid>
 </Domain>
</Xdmf>""" %(Nnodesb,Nnodesb, fn, nElementsb, nElementsb, fn, nElementsb, fn,nElementsb, fn)
   fid=open(fn+'.xdmf','w')
   fid.write(xdmf)
   fid.close()

def write_xdmfh5(fn, xyzb, tetrab, groupb, boundaryb):
   nNodesb=np.shape(xyzb)[0]
   nElementsb=np.shape(tetrab)[0]
   CreateXdmf(fn, Nnodesb, nElementsb)
   #Write h5 file
   import h5py
   h5f = h5py.File(fn,'w')
   h5f.create_dataset('/connect', tetrab.shape, dtype='i8')
   h5f['/connect'][:,:] = tetrab[:,:]
   h5f.create_dataset('/geometry', xyzb.shape, dtype='d')
   h5f['geometry'][:,:] = xyzb[:,:]
   h5f.create_dataset('/group', groupb.shape, dtype='i8')
   h5f['/group'][:] = groupb[:]
   h5f.create_dataset('/boundary', boundaryb.shape, dtype='i8')
   h5f['/boundary'][:] = boundaryb[:]
   h5f.close()
   print ("done writing %s.xdmf" %fn)


parser = argparse.ArgumentParser(description='Read hdf5 mesh and mirror it')
parser.add_argument('filename', help='filename of xdmf SeisSol mesh')
parser.add_argument('normal', nargs=1, help='0:x 1:y 2:z', type=int)
parser.add_argument('xc', nargs=1, help='coordinate on plane', type=float)
args = parser.parse_args()

xc = args.xc[0]
iN = args.normal[0]
sx = seissolxdmf.seissolxdmf(args.filename)
xyz = sx.ReadGeometry()
tetra = sx.ReadConnect()
nNodes=np.shape(xyz)[0]
nElements=np.shape(tetra)[0]

boundary = sx.ReadDataChunk('boundary', firstElement=0, nchunk=nElements, idt=0)
group = sx.ReadDataChunk('group', firstElement=0, nchunk=nElements, idt=0)

#Create Lookup and inverse Look table to relate old and new nodes
nodesLU = {}
inew = nNodes
for i in range(nNodes):
   if abs(xyz[i,iN]-xc)<1e-7:
     nodesLU[i]=i
   else:
     nodesLU[i]=inew
     inew = inew+1
invnodesLU = {v: k for k, v in nodesLU.items()}

Nnodesb = inew

#create new vertices
xyzb = np.zeros((inew,3))
xyzb[0:nNodes,:] = xyz
for inew in range(nNodes, Nnodesb):
   i = invnodesLU[inew]
   xyzb[inew,:] = xyzb[i,:] 
   xyzb[inew,iN] = xc-(xyzb[i,iN]-xc)
  

#create new connect
tetrab = np.zeros((nElements*2,4), dtype=int)
tetrab[0:nElements,:] = tetra
for i in range(nElements, 2*nElements):
   tetrab[i,:] = [nodesLU[k] for k in tetra[i-nElements,:]]
   #switch node 3 and 4 for getting a positive volume
   temp = tetrab[i,2]
   tetrab[i,2] = tetrab[i,3]
   tetrab[i,3] = temp

filename, file_extension = os.path.splitext(args.filename)
fn = filename+'_sym'

#create new group
groupb = np.zeros((nElements*2), dtype=int)
groupb[0:nElements] = group
for i in range(nElements, 2*nElements):
   groupb[i] = group[i-nElements]

boundaryb = np.zeros((nElements*2), dtype=int)
boundaryFace = np.zeros((nElements,4), dtype=int)

for faceId in range(0,4):
   boundaryFace[:,faceId] = (boundary >> (faceId*8)) & 0xFF;

boundaryFaceb = np.zeros((2*nElements,4))
boundaryFaceb[0:nElements,:]= boundaryFace

for ic in [1,3,5]:
   ids = np.where(boundaryFace==ic)
   elems=ids[0]
   faces=ids[1]
   #switch 0 and 1 for symetry
   ids0 = np.where(faces==0)
   ids1 = np.where(faces==1)
   faces[ids0]=1
   faces[ids1]=0
   for k in range(len(elems)):
      boundaryFaceb[nElements+elems[k],faces[k]]=ic

boundaryb = boundaryFaceb[:,0]+256*boundaryFaceb[:,1]+256*256*boundaryFaceb[:,2]+256*256*256*boundaryFaceb[:,3]

write_xdmfh5(fn, xyzb, tetrab, groupb, boundaryb)

