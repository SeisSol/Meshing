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
  [[nodiscard]] unsigned int numVertices() const { return puml.vertices().size(); };
  [[nodiscard]] unsigned int numElements() const { return puml.cells().size(); };
  [[nodiscard]] bool checkNeighbors() const;

  private:
  PUML::TETPUML puml;
  std::string_view bcToString(int id) const {
    switch(id) {
      case 0: return "regular";
      case 1: return "free surface";
      case 3: return "dynamic rupture";
      case 5: return "absorbing";
      default: return "";
    }
  }

};

#endif
