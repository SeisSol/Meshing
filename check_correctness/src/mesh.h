#ifndef MESH_H_
#define MESH_H_

#include <array>
#include <string>
#include <vector>

#include "PUML/PUML.h"

class Mesh {
  public:
  explicit Mesh(std::string const& fileName);
  ~Mesh() = default;
  [[nodiscard]] unsigned int numVertices() const { return puml.numOriginalVertices(); };
  [[nodiscard]] unsigned int numElements() const { return puml.numOriginalCells(); };
  [[nodiscard]] bool checkNeighbors() const;

  private:
  PUML::TETPUML puml;
};

#endif
