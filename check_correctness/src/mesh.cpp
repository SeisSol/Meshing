#include "mesh.h"

#include <mpi.h>

#include "PUML/PUML.h"
#include "PUML/Neighbor.h"
#include "PUML/Downward.h"

Mesh::Mesh(const std::string& fileName) {
  puml.open((fileName + ":/connect").c_str(), (fileName + ":/geometry").c_str());
  puml.addData((fileName + ":/group").c_str(), PUML::CELL);
  puml.addData((fileName + ":/boundary").c_str(), PUML::CELL);
  int numPartitionCells = puml.numOriginalCells();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Keep all cells on the same rank
  std::vector<int> partition;
  partition.reserve(numPartitionCells);
  for (size_t elemIdx = 0; elemIdx < numPartitionCells; elemIdx++) {
    partition[elemIdx] = rank;
  }
  puml.partition(partition.data());
  puml.generateMesh();
}

bool Mesh::checkNeighbors() const {
  bool result = true;
  unsigned int faceids[4]{};
  int neighbors[4]{};
  // check all cells
  for (size_t elemIdx = 0; elemIdx < numElements(); elemIdx++) {
    const auto& cell = puml.cells().at(elemIdx);
    const auto bndInfo = puml.cellData(1)[elemIdx];
    auto decodeBC = [&bndInfo](int side) { return (bndInfo >> (side*8)) & 0xFF; };

    PUML::Downward::faces(puml, cell, faceids);
    PUML::Neighbor::face(puml, elemIdx, neighbors);

    // check all sides
    for (size_t side = 0; side < 4; side++) {
      const auto& face = puml.faces()[faceids[side]];
      const auto sideBC = decodeBC(side);
      // if a face is a regular or a dr face, it has to have a neighbor on either this rank or somewhere else:
      if (sideBC == 0 || sideBC == 3) {
        if (neighbors[side] < 0 && !face.isShared()) {
          logInfo() << "Element" << cell.gid() << ", side" << side << " has a" << bcToString(sideBC) << "boundary condition, but the neighborig element doesn't exist";
          result = false;
        }
      }
      // absorbing or free surface boundaries must not have neighbor elements:
      else if (sideBC == 1 || sideBC == 5) {
        if (neighbors[side] >= 0 || face.isShared()) {
          logInfo() << "Element" << cell.gid() << ", side" << side << " has a" << bcToString(sideBC) << "boundary condition, but the neighborig element is not flagged -1";
          result = false;
        }
      }
      // ignore unknown boundary conditions and warn
      else {
        logWarning() << "Element" << cell.gid() << ", side" << side << " has a boundary condition, which I don't understand" << sideBC;
      }
    }
  }
  return result;
}
