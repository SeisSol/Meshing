// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Thomas Ulrich
// @author Sebastian Wolf

#include "Geometry.h"
#include <utils/logger.h>

void geometry::PointChecker::setupElements(const reader::Mesh& mesh) {
  elements.reserve(mesh.elementSize[0]);
  for (size_t elementIdx = 0; elementIdx < mesh.elementSize[0]; elementIdx++) {
    std::array<Eigen::Vector3d, 4> elementVertices;
    for (size_t j = 0; j < 4; j++) {
      auto vertexId = mesh.elementVertices[4 * elementIdx + j];
      elementVertices[j] = Eigen::Vector3d(&mesh.vertexCoordinates[3 * vertexId]);
    }
    elements.emplace_back(elementVertices);
  }
}

std::optional<int> geometry::PointChecker::pointInMesh(const Eigen::Vector3d& point) {
  std::optional<int> result = {};

#pragma omp parallel for
  for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
    if (elements[elementIdx].containsPoint(point)) {
      result = elementIdx;
    }
  }
  return result;
}

geometry::Element::Element(std::array<Eigen::Vector3d, 4> vertices)
    : vertices(std::move(vertices)) {
  for (size_t sideIdx = 0; sideIdx < 4; sideIdx++) {
    const auto& point0 = vertices[sideIdx];
    const auto& point1 = vertices[(sideIdx + 1) % 4];
    const auto& point2 = vertices[(sideIdx + 2) % 4];
    const auto& point3 = vertices[(sideIdx + 3) % 4];

    Eigen::Vector3d normal = (point1 - point0).cross(point2 - point0);
    normal.normalize();
    const auto dist = (point3 - point0).dot(normal);
    double const factor = -2 * static_cast<int>(dist > 0) + 1;
    faceNormals[sideIdx] = factor * normal;
  }
}

bool geometry::Element::containsPoint(const Eigen::Vector3d& point) {
  bool result = true;
  for (size_t sideIdx = 0; sideIdx < 4; sideIdx++) {
    const auto& point0 = vertices[sideIdx];
    const auto dist = (point - point0).dot(faceNormals[sideIdx]);
    result &= (dist < 0);
  }
  return result;
}
