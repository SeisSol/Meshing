// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Sebastian Wolf

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Reader.h"

#include <array>
#include <optional>
#include <vector>
#include <Eigen/Dense>

namespace geometry {
  class Element {
    public:
      Element(std::array<Eigen::Vector3d, 4> vertices);
      bool containsPoint(const Eigen::Vector3d& point);
    private:
      std::array<Eigen::Vector3d, 4> vertices;
      std::array<Eigen::Vector3d, 4> faceNormals;
  };

  class PointChecker {
    public:
      PointChecker(const reader::Mesh& mesh) {
        setupElements(mesh);
      }
      std::optional<int> pointInMesh(const Eigen::Vector3d& point);
    private:
      std::vector<Element> elements;
      void setupElements(const reader::Mesh& mesh);
  };
}


#endif
