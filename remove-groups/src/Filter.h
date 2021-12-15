#ifndef REMOVE_WATERLAYER_FILTER_H
#define REMOVE_WATERLAYER_FILTER_H

#include <cassert>
#include "Reader.h"

void removeGroups(Mesh &mesh,
                  const std::set<unsigned int> &groupsToRemove) {
  int removedCells = 0;
  // Map old element/vertex id -> new element/vertex id
  std::unordered_map<unsigned, unsigned> elementMap;
  std::unordered_map<unsigned, unsigned> vertexMap;

  std::cout << "Finding cells to remove" << std::endl;
  auto elementCounter = 0UL;
  auto vertexCounter = 0UL;
  for (auto cell = 0u; cell < mesh.elementSize; ++cell) {
    if (groupsToRemove.count(mesh.elementGroups[cell]) > 0) {
      // Remove cell
      ++removedCells;
      const auto &neighbors = mesh.neighbors[cell];
      // For all neighbors, set face type to free surface in slot of current cell
      for (auto neighborId : neighbors) {
        if (neighborId) {
          auto &neighborBnds = mesh.elementBoundaries[*neighborId];
          auto &neighborNeighbors = mesh.neighbors[*neighborId];

          // Find which nth neighbor we are of our neighbor
          auto faceId = -1;
          for (auto j = 0; j < 4; ++j) {
            auto &nn = neighborNeighbors[j];
            if (nn && *nn == cell) {
              faceId = j;
            }
          }
          neighborBnds[faceId] = 1; // Change to free surface bc

        }
      }
    } else {
      assert(elementMap.find(cell) == elementMap.end());
      elementMap[cell] = elementCounter;
      ++elementCounter;

      for (int i = 0; i < 4; ++i) {
        auto oldVertexId = mesh.connect[cell * 4 + i];
        if (vertexMap.find(oldVertexId) == vertexMap.end()) {
          vertexMap[oldVertexId] = vertexCounter++;
        }
      }

    }
    if (cell % (mesh.elementSize / 10) == 0) {
      std::cout << "Processed " << cell << "\t cells out of " << mesh.elementSize << std::endl;
    }
  }

  std::cout << "Going to remove " << removedCells << " cells." << std::endl;
  std::cout << "Begin removing" << std::endl;

  auto newElementGroups = std::vector<int>(elementCounter);
  auto newElementBoundaries = std::vector<Mesh::boundaries_t>(elementCounter);
  auto newElementConnect = std::vector<unsigned long>(elementCounter * 4);
  auto newElementNeighbors = std::vector<Mesh::neighbors_t>(elementCounter);

  for (auto cell = 0u; cell < mesh.elementSize; ++cell) {
    if (elementMap.count(cell) == 0) continue;
    auto newCellId = elementMap[cell];
    newElementBoundaries[newCellId] = mesh.elementBoundaries[cell];
    newElementGroups[newCellId] = mesh.elementGroups[cell];
    for (int i = 0; i < 4; ++i) {
      newElementConnect[newCellId * 4 + i] = vertexMap[mesh.connect[cell * 4 + i]];
    }
    for (int i = 0; i < 4; ++i) {
      auto oldNeighbor = mesh.neighbors[cell][i];
      if (oldNeighbor && elementMap.count(*oldNeighbor) > 0) {
        newElementNeighbors[newCellId][i] = elementMap[*oldNeighbor];
      } else {
        newElementNeighbors[newCellId][i] = std::nullopt;
      }
    }

    if (cell % (mesh.elementSize / 10) == 0) {
      std::cout << "Updated " << cell << "\t cells out of " << mesh.elementSize << std::endl;
    }
  }

  mesh.elementSize = elementCounter;

  mesh.elementBoundaries = std::move(newElementBoundaries);
  mesh.elementGroups = std::move(newElementGroups);
  mesh.connect = std::move(newElementConnect);
  mesh.neighbors = std::move(newElementNeighbors);

  std::cout << "Adjusting vertices" << std::endl;

  auto newVertices = std::vector<double>(vertexCounter * 3);
  for (unsigned vertex = 0; vertex < mesh.vertexSize; ++vertex) {
    if (vertexMap.find(vertex) == vertexMap.end()) continue;
    auto newVertexId = vertexMap[vertex];
    for (int i = 0; i < 3; ++i) {
      newVertices[newVertexId * 3 + i] = mesh.vertices[vertex * 3 + i];
    }
  }
  mesh.vertexSize = vertexCounter;
  mesh.vertices = std::move(newVertices);


}

#endif //REMOVE_WATERLAYER_FILTER_H
