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