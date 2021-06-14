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

#include "Reader.h"

#include <vector>
#include <optional>

#include "PUML/PUML.h"
#include "PUML/Neighbor.h"
#include "mpi.h"

Mesh::Mesh(const std::string &fileName) {
  PUML::TETPUML puml;
  puml.setComm(MPI_COMM_WORLD);
  puml.open((fileName + ":/connect").c_str(), (fileName + ":/geometry").c_str());
  puml.addData((fileName + ":/group").c_str(), PUML::CELL);
  puml.addData((fileName + ":/boundary").c_str(), PUML::CELL);
  puml.generateMesh();

  const auto &cells = puml.cells();
  const auto &faces = puml.faces();
  const auto &verticesPuml = puml.vertices();
  elementSize = cells.size();
  vertexSize = verticesPuml.size();

  elementBoundaries = std::vector<boundaries_t>(elementSize);
  elementGroups = std::vector<int>(elementSize);
  connect = std::vector<unsigned long>(elementSize * 4);
  auto connectOffset = 0;
  vertices = std::vector<double>(vertexSize * 3);
  auto verticesOffset = 0;

  neighbors = std::vector<neighbors_t>(elementSize);
  std::array<int, 4> neighborsRaw{};

  for (auto &vertex : verticesPuml) {
    const auto *coordinate = vertex.coordinate();
    for (int i = 0; i < 3; ++i) {
      vertices[verticesOffset++] = coordinate[i];
    }
  }
  for (int cell = 0; cell < elementSize; ++cell) {
    elementGroups[cell] = puml.cellData(0)[cell];

    std::array<unsigned int, 4> curVertices{};
    PUML::Downward::vertices(puml, cells[cell], curVertices.data());
    for (auto vertexId : curVertices) {
      connect[connectOffset++] = vertexId;
    }

    unsigned int faceids[4];
    PUML::Downward::faces(puml, cells[cell], faceids);

    PUML::Neighbor::face(puml, cell, neighborsRaw.data());
    neighbors_t curNeighbors{};
    for (int face = 0; face < 4; ++face) {
      curNeighbors[face] = (neighborsRaw[face] >= 0) ? std::optional(neighborsRaw[face]) : std::nullopt;
      elementBoundaries[cell][face] = decodeBoundaryCondition(puml.cellData(1)[cell], face);
    }
    neighbors[cell] = curNeighbors;
  }


}

int decodeBoundaryCondition(int encoded, int faceId) {
  return (encoded >> (faceId * 8)) & 0xFF;
}