#include "mesh.h"

#include "PUML/PUML.h"

Mesh::Mesh(const std::string& fileName) {
  PUML::TETPUML puml;
  puml.open((fileName + ":/connect").c_str(), (fileName + ":/geometry").c_str());
  puml.addData((fileName + ":/group").c_str(), PUML::CELL);

  const unsigned int nElements = puml.numOriginalCells();
  const unsigned int nVertex = puml.numOriginalVertices();
  elementVertices.reserve(nElements);
  vertexCoordinates.reserve(nVertex);
  elementGroups.reserve(nElements);

  for (unsigned int i = 0; i < nVertex; i++) {
    const auto vertexCoordinate = puml.originalVertices()[i];
    std::array<double, 3> tmp({vertexCoordinate[0], vertexCoordinate[1], vertexCoordinate[2]});
    vertexCoordinates.emplace_back(tmp);
  }

  for (unsigned int i = 0; i < nElements; i++) {
    const auto elementVertex = puml.originalCells()[i];
    std::array<unsigned long, 4> tmp({elementVertex[0], elementVertex[1], elementVertex[2], elementVertex[3]});
    elementVertices.emplace_back(tmp);
  }

  for (unsigned int i = 0; i < nElements; i++) {
    elementGroups.emplace_back(puml.cellData(0)[i]);
  }
}

std::vector<std::array<double, 3>> Mesh::getElementBarycenters() const {
  std::vector<std::array<double, 3>> barycenters;
  constexpr auto numOfVertices = 4.0;
  barycenters.reserve(elementVertices.size());
  for (const auto& vertexIndices : elementVertices) {
    std::array<double, 3> b{};
    for (const auto& vertexId : vertexIndices) {
      for (size_t dim = 0; dim < 3; dim++) {
        b.at(dim) += vertexCoordinates[vertexId][dim] / numOfVertices;
      }
    }
    barycenters.emplace_back(b);
  }

  return barycenters;
}
