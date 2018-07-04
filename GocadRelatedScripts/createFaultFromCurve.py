import numpy as np
from math import tan,pi
from scipy.interpolate import interp1d
import os

import argparse
parser = argparse.ArgumentParser(description='create curved fault geometry from pl file')
parser.add_argument('filename', help='fault trace (*.pl) or ascii file (2 or 3 columns)')
parser.add_argument('dip', help='dip value or name of ascii file with 2 columns depth dip')
parser.add_argument('--constV1',  dest='constV1', action='store_true' , default = False, help='dip horizontal direction from 2 extreme points of pl file (else changing along trace)')
parser.add_argument('--translate', nargs=2, metavar=('x0', 'y0'), default = ([0,0]), help='translates all nodes by (x0,y0)', type=float)
parser.add_argument('--dd', nargs=1, metavar=('dd'), default = ([1e3]), help='sampling along depth', type=float)
parser.add_argument('--maxdepth', nargs=1, metavar=('maxdepth'), default = ([20e3]), help='max depth (positive) of fault', type=float)
parser.add_argument('--extend', nargs=1, metavar=('extend'), default = ([00e3]), help='extend toward z= extend (positive)', type=float)
args = parser.parse_args()



constantdip=True
dx = args.dd[0]

try:
   dip = float(args.dip)*pi/180.
except ValueError:
   constantdip=False
   print(args.dip)
   depthdip=np.loadtxt(args.dip)
   deptha = depthdip[:,0]
   dipa = depthdip[:,1]*pi/180.
   dipangle = interp1d(deptha, dipa, kind='linear')

bn = os.path.basename(args.filename)
ext = bn.split('.')[1]
if ext=='pl':
   nodes=[]
   fid = open(args.filename)
   lines = fid.readlines()
   fid.close()
   for li in lines:
      if li.startswith('VRTX'):
          lli = li.split()
          nodes.append([float(lli[2]), float(lli[3]), float(lli[4])])

   nodes= np.asarray(nodes)
else:
   nodes=np.loadtxt(args.filename)
   ndim = np.shape(nodes)[1]
   if ndim==2:
      nx = np.shape(nodes)[0]
      b = np.zeros((nx,1))
      nodes = np.append(nodes, b, axis=1)

nodes[:,0] = nodes[:,0]+args.translate[0]
nodes[:,1] = nodes[:,1]+args.translate[1]

print(nodes)
maxdepth = -args.maxdepth[0]
depth = np.arange(maxdepth,0.,dx)

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
         mydip = dipangle(depth[nd-j])
      vertices[i,j,:] = vertices[i,j-1,:] + dx * (1./tan(mydip)* v1 -uz)

if args.extend[0]>0:
   #Extension toward z plus
   depth = np.arange(0,args.extend[0],dx)
   print(depth)
   nd2 = np.shape(depth)[0]
   vertices2 = np.zeros((nx,nd2,3))
   vertices2[:,0,:] = nodes
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
      for j in range(1,nd2):
         if constantdip:
            mydip = dip
         else:
            mydip = dipangle(depth[nd2-j])
         vertices2[i,j,:] = vertices2[i,j-1,:] - dx * (1./tan(mydip)* v1 -uz)

bn = os.path.basename(args.filename)
prefix = bn.split('.')[0]
fid=open(prefix + '0.dat','w')
for j in reversed(range(0,nd)):
   for i in range(nx):
       fid.write('%15.10e %15.10e %15.10e\n' %(vertices[i,j,0],vertices[i,j,1],vertices[i,j,2]))

if args.extend[0]>0:
   for j in range(1,nd2):
      for i in range(nx):
          fid.write('%15.10e %15.10e %15.10e\n' %(vertices2[i,j,0],vertices2[i,j,1],vertices2[i,j,2]))
fid.close()

print('done! NX=%d' %np.shape(vertices)[0])



