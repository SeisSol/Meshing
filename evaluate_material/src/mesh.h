#ifndef MESH_H_
#define MESH_H_

#include <array>
#include <string>
#include <vector>

class Mesh {
  public:
  explicit Mesh(std::string const& fileName);
  ~Mesh() = default;
  [[nodiscard]] std::vector<std::array<double, 3>> getElementBarycenters() const;
  [[nodiscard]] const std::vector<std::array<double, 3>>& getVertexCoordinates() const { return vertexCoordinates; };
  [[nodiscard]] const std::vector<std::array<unsigned long, 4>>& getElementVertices() const { return elementVertices; };
  [[nodiscard]] const std::vector<int>& getElementGroups() const { return elementGroups; };

  private:
  std::vector<std::array<unsigned long, 4>> elementVertices;
  std::vector<std::array<double, 3>> vertexCoordinates;
  std::vector<int> elementGroups;
};

#endif
