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
#include "LinearReceiverSearch.h"

#include <array>
#include <cmath>
#include <limits>

#ifdef USE_NETCDF
int const FACE2NODES[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
#else
int const FACE2NODES[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
#endif

void setElevation(int partition, double depth, Mesh const& mesh,
                  LinearReceiverSearch & receiver) {
  for (unsigned element = 0; element < mesh.elementSize[partition]; ++element) {
    for (unsigned face = 0; face < 4; ++face) {
      // Consider only free surface boundaries
      if (mesh.elementBoundaries[4*element + face] == 1) {
        std::array<Point, 3> vertices{};
        for (unsigned node = 0; node < 3; ++node) {
          double* vbegin = &mesh.vertexCoordinates[3*mesh.elementVertices[ 4*element + FACE2NODES[face][node] ]];
          vertices[node] = {
             *(vbegin),
             *(vbegin + 1),
             *(vbegin + 2),
          };
        }
        Triangle triangle{vertices[0], vertices[1], vertices[2]};
        receiver.applyFunctionOnAllReceiversInTriangle(triangle, [&triangle, depth](Point& receiver)
        {
          double local_depth = depth;
          if (!std::isnan(receiver.z)) {
            local_depth = std::max(depth, std::fabs(receiver.z));
          }
          receiver.z = (triangle.faceDist - triangle.faceNormal[0] * receiver.x - triangle.faceNormal[1] * receiver.y) / triangle.faceNormal[2];
          if (triangle.faceNormal[2] >= 0) {
            receiver.z -= local_depth;
          } else {
            receiver.z += local_depth;
          }
          return receiver;
        });
      }
    }
  }
}
