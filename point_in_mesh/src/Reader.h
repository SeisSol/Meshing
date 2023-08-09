// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Sebastian Wolf

#ifndef READER_H_
#define READER_H_

#include <string>
#include <vector>
#include <Eigen/Dense>

namespace reader {
  class Mesh {
    public:
      explicit Mesh(std::string const& fileName);
      int partitions;
      std::vector<int> elementSize;
      std::vector<int> vertexSize;
      std::vector<int> elementBoundaries;
      std::vector<int> elementVertices;
      std::vector<double> vertexCoordinates;
  };

  std::vector<Eigen::Vector3d> readReceiverFile(std::string const& fileName);
}

#endif
