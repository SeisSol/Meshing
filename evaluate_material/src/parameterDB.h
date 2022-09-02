#ifndef PARAMETER_DB_H_
#define PARAMETER_DB_H_

#include <string>

#include "easi/Component.h"
#include "easi/Query.h"

class ParameterDB {
  public:
  ParameterDB(const std::string& fileName);
  ~ParameterDB();
  [[nodiscard]] easi::Query generateQuery(const std::vector<std::array<double, 3>>& points,
                                          const std::vector<int>& groups) const;
  std::pair<std::vector<std::string>, std::vector<std::vector<double>>> evaluate(easi::Query& query) const;

  private:
  easi::Component* model;
  static constexpr unsigned int DIM = 3;
};

#endif
