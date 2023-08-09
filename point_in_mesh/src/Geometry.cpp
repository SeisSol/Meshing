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
      auto vertexId = mesh.elementVertices[4*elementIdx+j];
      elementVertices[j] = Eigen::Vector3d(&mesh.vertexCoordinates[3*vertexId]);
    }
    const auto center = (elementVertices[0] + elementVertices[1] + elementVertices[2] + elementVertices[3]) * 0.25; 
    elements.emplace_back(Element(elementVertices));
  }
}

std::optional<int> geometry::PointChecker::pointInMesh(const Eigen::Vector3d& point) {
  for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
    if (elements[elementIdx].containsPoint(point)) {
      return elementIdx;
    }
  }
  return {};
}

geometry::Element::Element(std::array<Eigen::Vector3d, 4> vertices) : vertices(std::move(vertices)) {
  for (size_t sideIdx = 0; sideIdx < 4; sideIdx++) {
    const auto& p0 = vertices[sideIdx];
    const auto& p1 = vertices[(sideIdx+1)%4];
    const auto& p2 = vertices[(sideIdx+2)%4];
    const auto& p3 = vertices[(sideIdx+3)%4];

    Eigen::Vector3d normal = (p1 - p0).cross(p2-p0);
    normal.normalize();
    const auto dist = (p3 - p0).dot(normal);
    double factor = -2 * (dist > 0) + 1;
    faceNormals[sideIdx] = factor * normal;
  }
}

bool geometry::Element::containsPoint(const Eigen::Vector3d& point) {
  bool result = true;
  for (size_t sideIdx = 0; sideIdx < 4; sideIdx++) {
    const auto& p0 = vertices[sideIdx];
    const auto dist = (point - p0).dot(faceNormals[sideIdx]);
    result &= (dist < 0);
  }
  return result;
}
