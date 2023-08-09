// SPDX-License-Identifier: BSD-3-Clause
// @author Carsten Uphoff
// @author Thomas Ulrich 
// @author Sebastian Wolf

#include "Reader.h"

#include <utils/logger.h>
#include <PUML/PUML.h>

//#include <cstdlib>
//#include <cmath>
//#include <iostream>
//#include <iomanip>
#include <fstream>
//#include <sstream>

std::vector<Eigen::Vector3d> reader::readReceiverFile(std::string const& fileName) {
  using Vector3d = Eigen::Vector3d;
  std::vector<Vector3d> locations;
  std::ifstream in(fileName.c_str());
  std::string line;

  while (std::getline(in, line)) {
    std::istringstream iss(line);
    Vector3d p;
    int coord = 0;
    while (coord < 3 && iss.good()) {
      iss >> p(coord++);
    }
    locations.push_back(p);
    if (iss.bad()) {
      logError() << "An error occurred while reading the receiver file.";
    }
  }
  
  return locations;
}

//void writeReceiverFile(KDTree const& tree, std::string const& fileName) {
//  std::ofstream out(fileName.c_str());
//  out << std::scientific << std::setprecision(16);
//  
//  int failureCounter = 0;
//  Point const* points = tree.points();
//  std::vector<Point> sortedPoints(tree.numPoints());
//  for (unsigned p = 0; p < tree.numPoints(); ++p) {
//    sortedPoints[tree.index(p)] = points[p];
//  }
//  for (auto const& point : sortedPoints) {
//    if (!std::isnan(point.z)) {
//      out << point.x << " " << point.y << " " << point.z << std::endl;
//    } else {
//      std::cerr << "Warning: Did not find elevation for receiver at (" << point.x << ", " << point.y << ")." << std::endl;
//      ++failureCounter;
//    }
//  }
//
//  if (failureCounter > 0) {
//    std::cerr << failureCounter << " points were not written due to missing elevation." << std::endl;
//  }
//}
//
//#ifdef USE_NETCDF
//void check_err(const int stat, const int line, const char *file) {
//  if (stat != NC_NOERR) {
//    fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
//    fflush(stderr);
//    exit(-1);
//  }
//}
//#endif

reader::Mesh::Mesh(std::string const& fileName){
  PUML::TETPUML puml;
  puml.open((fileName + ":/connect").c_str(), (fileName + ":/geometry").c_str());
  puml.addData((fileName + ":/boundary").c_str(), PUML::CELL);

  auto nElements = puml.numOriginalCells();
  auto nVertex = puml.numOriginalVertices();
  elementVertices = std::vector<int>(nElements*4);
  vertexCoordinates = std::vector<double>(nVertex*3);
  elementBoundaries = std::vector<int>(nElements*4);

  partitions = 1;
  elementSize = std::vector<int>(partitions);
  vertexSize = std::vector<int>(partitions);
  elementSize[0] = nElements;
  vertexSize[0]= nVertex;

  double* vertexCoordinate = new double[3]; 
  for (unsigned int i = 0; i < nVertex; i++) {
    vertexCoordinate = (double*) puml.originalVertices()[i];
    for (unsigned int j = 0; j < 3; j++) {
       vertexCoordinates[3*i+j] =  vertexCoordinate[j];
    }
  }

  unsigned long* elementVertex = new unsigned long[4]; 
  for (unsigned int i = 0; i < nElements; i++) {
    elementVertex = (unsigned long *) puml.originalCells()[i];
    for (unsigned int j = 0; j < 4; j++) {
       elementVertices[4*i+j] =  (int) elementVertex[j];
    }
  }

  const int* boundary = puml.cellData(0);
  for (unsigned int i = 0; i < nElements; i++) {
    for (unsigned int face = 0; face < 4; face++) {
       elementBoundaries[4*i + face] = (boundary[i] >> (face*8)) & 0xFF;
    }
  }
}
