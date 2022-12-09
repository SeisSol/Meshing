#include "mesh.h"

#include <numeric>

#include <mpi.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include "PUML/Neighbor.h"
#include "PUML/Downward.h"
#pragma clang diagnostic pop

Mesh::Mesh(const std::string& fileName) {
  puml.open((fileName + ":/connect").c_str(), (fileName + ":/geometry").c_str());
  puml.addData((fileName + ":/group").c_str(), PUML::CELL);
  puml.addData((fileName + ":/boundary").c_str(), PUML::CELL);
  const auto numTotalCells = puml.numTotalCells();
  std::vector<int> cellIdsAsInFile;
  cellIdsAsInFile.reserve(numTotalCells);
  std::iota(cellIdsAsInFile.begin(), cellIdsAsInFile.begin(), 0);
  puml.addData(cellIdsAsInFile.data(), numTotalCells, PUML::CELL);
  
  // Keep all cells on the same rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::size_t numPartitionCells = puml.numOriginalCells();
  std::vector<int> partition;
  partition.reserve(numPartitionCells);
  for (std::size_t elemIdx = 0; elemIdx < numPartitionCells; elemIdx++) {
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
  for (size_t elemIdx = 0; elemIdx < this->numElements(); elemIdx++) {
    const auto& cell = puml.cells().at(elemIdx);
    const auto bndInfo = puml.cellData(1)[elemIdx];
    const auto cellIdAsInFile = puml.cellData(2)[elemIdx];
    auto decodeBC = [&bndInfo](std::size_t side) { return (bndInfo >> (side*8)) & 0xFF; };

    PUML::Downward::faces(puml, cell, faceids);
    PUML::Neighbor::face(puml, elemIdx, neighbors);

    // check all sides
    for (size_t side = 0; side < 4; side++) {
      const auto& face = puml.faces()[faceids[side]];
      const auto sideBC = decodeBC(side);
      // if a face is an internal face, it has to have a neighbor on either this rank or somewhere else:
      if (bcToType(sideBC) == BCType::internal) {
        if (neighbors[side] < 0 && !face.isShared()) {
          logInfo() << "Element" << cellIdAsInFile << ", side" << side << " has a" << bcToString(sideBC) << "boundary condition, but the neighboring element doesn't exist";
          result = false;
        }
      }
      // external boundaries must not have neighboring elements:
      else if (bcToType(sideBC) == BCType::external) {
        if (neighbors[side] >= 0 || face.isShared()) {
          logInfo() << "Element" << cellIdAsInFile << ", side" << side << " has a" << bcToString(sideBC) << "boundary condition, but the neighboring element is not flagged -1";
          result = false;
        }
      }
      // ignore unknown boundary conditions and warn
      else {
        logWarning() << "Element" << cellIdAsInFile << ", side" << side << " has a boundary condition, which I don't understand" << sideBC;
      }
    }
  }
  return result;
}
