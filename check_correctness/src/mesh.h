#ifndef MESH_H_
#define MESH_H_

#include <array>
#include <map>
#include <string>
#include <vector>

#include "PUML/PUML.h"


class Mesh {
  public:
  explicit Mesh(std::string const& fileName);
  ~Mesh() = default;
  [[nodiscard]] std::size_t numVertices() const { return puml.vertices().size(); }
  [[nodiscard]] std::size_t numElements() const { return puml.cells().size(); }
  [[nodiscard]] bool checkNeighbors() const;

  private:
  PUML::TETPUML puml;
  std::string bcToString(int id) const {
    switch(id) {
      case 0: return std::string("regular");
      case 1: return std::string("free surface");
      case 3: return std::string("dynamic rupture");
      case 5: return std::string("absorbing");
      default: return std::string("");
    }
  }

};

#endif
