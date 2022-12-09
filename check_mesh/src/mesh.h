#ifndef MESH_H_
#define MESH_H_

#include <array>
#include <map>
#include <sstream>
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

  enum class BCType {
    internal, external, unknown
  };

  BCType bcToType(int id) const {
    if (id == 0 || id == 3 || id > 64) {
      return BCType::internal;
    } else if (id == 1 || id == 5 || id == 6 || id == 6) {
      return BCType::external;
    } else {
      return BCType::unknown;
    }
  }

  std::string bcToString(int id) const {
    if (id == 0) { return std::string("regular"); }
    else if (id == 1) { return std::string("free surface"); }
    else if (id == 3) { return std::string("dynamic rupture"); }
    else if (id == 5) { return std::string("absorbing"); }
    else if (id == 6) { return std::string("periodic"); }
    else if (id > 64) { 
      std::stringstream s;
      s << "fault-tagging (" << id << ")";
      return s.str();
    } else {return std::string(""); }
  }

};

#endif
