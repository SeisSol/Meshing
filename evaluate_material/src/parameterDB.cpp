#include "parameterDB.h"

#include <easi/Component.h>
#include <easi/Query.h>
#include <easi/ResultAdapter.h>
#include <easi/YAMLParser.h>

ParameterDB::ParameterDB(const std::string& fileName) {
  auto parser = easi::YAMLParser(DIM);
  model = parser.parse(fileName);
}

ParameterDB::~ParameterDB() { delete model; }

easi::Query ParameterDB::generateQuery(const std::vector<std::array<double, DIM>>& points,
                                       const std::vector<int>& groups) const {
  unsigned nPoints = points.size();
  auto query = easi::Query(nPoints, DIM);
  for (unsigned pointId = 0; pointId < nPoints; pointId++) {
    query.group(pointId) = groups[pointId];
    for (unsigned dim = 0; dim < DIM; dim++) {
      query.x(pointId, dim) = points[pointId][dim];
    }
  }
  return query;
}

std::pair<std::vector<std::string>, std::vector<std::vector<double>>> ParameterDB::evaluate(easi::Query& query) const {
  auto supplied = model->suppliedParameters();
  auto parameters = std::vector<std::string>(supplied.begin(), supplied.end());
  auto adapter = easi::ArraysAdapter<>{};

  auto material = std::vector<std::vector<double>>(parameters.size());
  auto it = material.begin();
  for (auto const& p : parameters) {
    it->resize(query.numPoints());
    adapter.addBindingPoint(p, it->data());
    ++it;
  }
  model->evaluate(query, adapter);
  return {parameters, material};
}
