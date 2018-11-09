/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Thomas Ulrich 
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "Geometry.h"

#include <limits>

#ifdef USE_NETCDF
int const FACE2NODES[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
#else
int const FACE2NODES[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
#endif

struct Support {
  // limits[0][:] = (min x, max x)
  // limits[1][:] = (min y, max y)
  double limits[2][2];

  Support() {
    limits[0][0] = std::numeric_limits<double>::max();
    limits[0][1] = std::numeric_limits<double>::min();
    limits[1][0] = std::numeric_limits<double>::max();
    limits[1][1] = std::numeric_limits<double>::min();
  }
  
  double operator()(int splitdim, int side) const {
    return limits[splitdim][side];
  }
};

struct Action {
  Point vertices[3];
  double normals[3][2];
  double dist[3];
  double faceNormal[3];
  double faceDist;
  double depth;
  
  void determineNormals() {
    for (unsigned side = 0; side < 3; ++side) {
      double* v0 = vertices[side].coords;
      double* v1 = vertices[(side+1)%3].coords;
      double* vtest = vertices[(side+2)%3].coords;
      // normal = (y,-x)
      normals[side][0] = v1[1]-v0[1];
      normals[side][1] = -v1[0]+v0[0];

      dist[side] = normals[side][0] * v0[0] + normals[side][1] * v0[1];
      
      if (normals[side][0] * vtest[0] + normals[side][1] * vtest[1] > dist[side]) {
        normals[side][0] *= -1.;
        normals[side][1] *= -1.;
      }
    }
    
    double a[3];
    double b[3];
    for (unsigned d = 0; d < 3; ++d) {
      a[d] = vertices[1].coords[d]-vertices[0].coords[d];
      b[d] = vertices[2].coords[d]-vertices[0].coords[d];
    }
    faceNormal[0] = a[1]*b[2] - a[2]*b[1];
    faceNormal[1] = a[2]*b[0] - a[0]*b[2];
    faceNormal[2] = a[0]*b[1] - a[1]*b[0];
    faceDist = faceNormal[0] * vertices[0].coords[0] + faceNormal[1] * vertices[0].coords[1] + faceNormal[2] * vertices[0].coords[2];
  }
  
  double operator()(Point& receiver) {
    bool inside = true;
    for (unsigned side = 0; side < 3; ++side) {
      inside = inside && (normals[side][0] * receiver.x + normals[side][1] * receiver.y <= dist[side]);
    }
    
    if (inside) {
      receiver.z = (faceDist - faceNormal[0] * receiver.x - faceNormal[1] * receiver.y) / faceNormal[2];
      if (faceNormal[2] >= 0) {
        receiver.z -= depth;
      } else {
        receiver.z += depth;
      }
    }
  }
};

void setElevation(int partition, double depth, Mesh const& mesh, KDTree& tree) {
  for (unsigned element = 0; element < mesh.elementSize[partition]; ++element) {
    for (unsigned face = 0; face < 4; ++face) {
      // Consider only free surface boundaries
      if (mesh.elementBoundaries[4*element + face] == 1) {
        Action act;
        act.depth = depth;
        Support sup;
        for (unsigned node = 0; node < 3; ++node) {
          double* vbegin = &mesh.vertexCoordinates[3*mesh.elementVertices[ 4*element + FACE2NODES[face][node] ]];
          act.vertices[node].coords[0] = *(vbegin);
          act.vertices[node].coords[1] = *(vbegin + 1);
          act.vertices[node].coords[2] = *(vbegin + 2);
          sup.limits[0][0] = std::min(sup.limits[0][0], act.vertices[node].coords[0]);
          sup.limits[0][1] = std::max(sup.limits[0][1], act.vertices[node].coords[0]);
          sup.limits[1][0] = std::min(sup.limits[1][0], act.vertices[node].coords[1]);
          sup.limits[1][1] = std::max(sup.limits[1][1], act.vertices[node].coords[1]);
        }
        act.determineNormals();
        tree.search(sup, act);
      }
    }
  }
}
