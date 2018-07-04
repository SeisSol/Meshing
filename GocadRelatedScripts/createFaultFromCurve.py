import numpy as np
from math import tan,pi
from scipy.interpolate import interp1d
import configparser

import argparse
parser = argparse.ArgumentParser(description='create curved fault geometry from pl file')
parser.add_argument('filename', help='fault trace (*.pl)')
parser.add_argument('dip', help='dip value or name of ascii file with 2 columns depth dip')
parser.add_argument('--constV1',  dest='constV1', action='store_true' , default = False, help='dip horizontal direction from 2 extreme points of pl file (else changing along trace)')
args = parser.parse_args()



constantdip=True

try:
   dip = float(args.dip)*pi/180.
except ValueError:
   constantdip=False
   print(args.dip)
   depthdip=np.loadtxt(args.dip)
   deptha = depthdip[:,0]
   dipa = depthdip[:,1]*pi/180.
   dipangle = interp1d(deptha, dipa, kind='linear')

nodes=[]
fid = open(args.filename)
lines = fid.readlines()
fid.close()
for li in lines:
   if li.startswith('VRTX'):
       lli = li.split()
       nodes.append([float(lli[2]), float(lli[3]), float(lli[4])])

nodes= np.asarray(nodes)
print(nodes)
dx=1e3
maxdepth=-25e3
depth = np.arange(maxdepth,0, dx)

nx = np.shape(nodes)[0]
nd = np.shape(depth)[0]
vertices = np.zeros((nx,nd,3))
uz = np.array([0, 0, 1])
vertices[:,0,:] = nodes

for i in range(0,nx):
   if not args.constV1:
      if i+1!=nx:
         v0 = vertices[i+1,0,:]-vertices[i,0,:]
      else:
         v0 = vertices[i,0,:]-vertices[i-1,0,:]
   else:
      v0 = vertices[nx-1,0,:]-vertices[0,0,:]

   v0[2] = 0
   v0 = v0/np.linalg.norm(v0)
   v1 = np.array([-v0[1], v0[0], 0])      
   for j in range(1,nd):
      if constantdip:
         mydip = dip
      else:
         mydip = dipangle(depth[j])
      vertices[i,j,:] = vertices[i,j-1,:] + dx * (tan(mydip)* v1 -uz)

import os
bn = os.path.basename(args.filename)
prefix = bn.split('.')[0]
fid=open(prefix + '.dat','w')
for i in range(nx):
   for j in range(0,nd):
       fid.write('%15.10e %15.10e %15.10e\n' %(vertices[i,j,0],vertices[i,j,1],vertices[i,j,2]))
fid.close()



print(np.shape(vertices))
