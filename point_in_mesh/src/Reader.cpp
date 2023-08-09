// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Thomas Ulrich
// @author Sebastian Wolf

#include "Reader.h"

#include <PUML/PUML.h>
#include <utils/logger.h>

#include <fstream>

/*
 * Reads a SeisSol receiver file (assumes 3 floating point values per row).
 * Returns a Vector of the points.
 * @param fileName receiver file name
 * @return vector containing all receivers in the same order as in the file.
 */
std::vector<Eigen::Vector3d> reader::readReceiverFile(std::string const& fileName) {
  using Vector3d = Eigen::Vector3d;
  std::vector<Vector3d> locations;
  std::ifstream inStream(fileName.c_str());
  std::string line;

  while (std::getline(inStream, line)) {
    std::istringstream iss(line);
    Vector3d point;
    int coord = 0;
    while (coord < 3 && iss.good()) {
      iss >> point(coord++);
    }
    locations.push_back(point);
    if (iss.bad()) {
      logError() << "An error occurred while reading the receiver file.";
    }
  }

  return locations;
}

reader::Mesh::Mesh(std::string const& fileName) {
  PUML::TETPUML puml;
  puml.open((fileName + ":/connect").c_str(), (fileName + ":/geometry").c_str());
  puml.addData((fileName + ":/boundary").c_str(), PUML::CELL);

  auto nElements = puml.numOriginalCells();
  auto nVertex = puml.numOriginalVertices();
  elementVertices = std::vector<int>(nElements * 4);
  vertexCoordinates = std::vector<double>(nVertex * 3);
  elementBoundaries = std::vector<int>(nElements * 4);

  partitions = 1;
  elementSize = std::vector<int>(partitions);
  vertexSize = std::vector<int>(partitions);
  elementSize[0] = nElements;
  vertexSize[0] = nVertex;

  auto* vertexCoordinate = new double[3];
  for (unsigned int i = 0; i < nVertex; i++) {
    vertexCoordinate = (double*)puml.originalVertices()[i];
    for (unsigned int j = 0; j < 3; j++) {
      vertexCoordinates[3 * i + j] = vertexCoordinate[j];
    }
  }

  auto* elementVertex = new unsigned long[4];
  for (unsigned int i = 0; i < nElements; i++) {
    elementVertex = (unsigned long*)puml.originalCells()[i];
    for (unsigned int j = 0; j < 4; j++) {
      elementVertices[4 * i + j] = (int)elementVertex[j];
    }
  }

  const int* boundary = puml.cellData(0);
  for (unsigned int i = 0; i < nElements; i++) {
    for (unsigned int face = 0; face < 4; face++) {
      elementBoundaries[4 * i + face] = (boundary[i] >> (face * 8)) & 0xFF;
    }
  }
}
