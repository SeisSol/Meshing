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
  //logInfo() << "Keep all" << numPartitionCells << "elements on rank" << rank;
  std::vector<int> partition;
  partition.reserve(numPartitionCells);
  for (size_t elemIdx = 0; elemIdx < numPartitionCells; elemIdx++) {
    partition[elemIdx] = rank;
  }
  puml.partition(partition.data());
  //logInfo() << "Done.";
  //logInfo() << "Generate mesh";
  puml.generateMesh();
  //logInfo() << "Done.";
}

bool Mesh::checkNeighbors() const {
  const size_t numElements = this->numElements();
  for (size_t elemIdx = 0; elemIdx < numElements; elemIdx++) {
    const auto& cell = puml.cells().at(elemIdx);

    const auto bndInfo = puml.cellData(1)[elemIdx];
    unsigned int faceids[4];
    PUML::Downward::faces(puml, cell, faceids);
    int neighbors[4];
    PUML::Neighbor::face(puml, elemIdx, neighbors);
    auto decodeBC = [&bndInfo](int side) { return (bndInfo >> (side*8)) & 0xFF; };

    //logInfo() << "Element" << elemIdx;
    for (size_t side = 0; side < 4; side++) {
      const auto& face = puml.faces()[faceids[side]];
      //logInfo() << "side" << side << ":" << decodeBC(side) << ", " << neighbors[side] << face.isShared();
      // if a face is a regular or a dr face, it has to have a neighbor on either this rank or somewhere else:
      if (decodeBC(side) == 0 || decodeBC(side) == 3) {
        if (neighbors[side] < 0 && !face.isShared()) {
          logInfo() << "Element" << cell.gid() << ", side" << side << " has a regular or DR boundary condition, but the neighborig element doesn't exist";
          return false;
        }
      }
      // absorbing or free surface boundaries must not have neighbor elements:
      else {
        if (neighbors[side] >= 0 || face.isShared()) {
          logInfo() << "Element" << cell.gid() << ", side" << side << " has an absorbing or free surface boundary conditin, but the neighborig element is not flagged -1";
          return false;
        }
      }
    }
  }
  return true;
}
