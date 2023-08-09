// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Sebastian Wolf

#ifndef READER_H_
#define READER_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace reader {
/*
 * Wraps a PUML mesh
 */
class Mesh {
  public:
  /*
   * Initializes a PUML mesh from an hdf5 file.
   * @param fileName mesh file name
   */
  explicit Mesh(std::string const& fileName);
  int partitions;
  std::vector<int> elementSize;
  std::vector<int> vertexSize;
  std::vector<int> elementBoundaries;
  std::vector<int> elementVertices;
  std::vector<double> vertexCoordinates;
};

/*
 * Reads a SeisSol receiver file (assumes 3 floating point values per row).
 * Returns a Vector of the points.
 * @param fileName receiver file name
 * @return vector containing all receivers in the same order as in the file.
 */
std::vector<Eigen::Vector3d> readReceiverFile(std::string const& fileName);
} // namespace reader

#endif
