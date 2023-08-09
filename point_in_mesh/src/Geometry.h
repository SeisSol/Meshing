// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Sebastian Wolf

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Reader.h"

#include <Eigen/Dense>
#include <array>
#include <optional>
#include <vector>

namespace geometry {
/*
 * Represents a tetrahedral element
 */
class Element {
  public:
  /*
   * Initializes an element from its four vertices.
   * In particular, it computes the four outward pointing normals.
   * @param vertices Array of the four vertices
   */
  Element(std::array<Eigen::Vector3d, 4> vertices);
  /*
   * Checks whether a point lies within the element.
   * @param point point to check.
   */
  bool containsPoint(const Eigen::Vector3d& point);

  private:
  std::array<Eigen::Vector3d, 4> vertices;
  std::array<Eigen::Vector3d, 4> faceNormals;
};

class PointChecker {
  public:
  /*
   * Creates a point checker based on a mesh
   * @param mesh reference to the mesh
   */
  PointChecker(const reader::Mesh& mesh) { setupElements(mesh); }
  /*
   * Checks whethe a point lies within the mesh.
   * @param point point to check.
   * @return if the point lies in the mesh, return the element id. If the
   * point does not lie in the mesh, return an empty optional.
   */
  std::optional<int> pointInMesh(const Eigen::Vector3d& point);

  private:
  std::vector<Element> elements;
  void setupElements(const reader::Mesh& mesh);
};
} // namespace geometry

#endif
