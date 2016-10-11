/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "Reader.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <netcdf.h>

std::vector<Point> readReceiverFile(std::string const& fileName) {
  std::vector<Point> locations;
  std::ifstream in(fileName.c_str());
  std::string line;

  while (std::getline(in, line)) {
    std::istringstream iss(line);
    Point p;
    iss >> p.x >> p.y;
    p.z = NAN;
    locations.push_back(p);
    if (iss.bad()) {
      std::cerr << "An error occurred while reading the receiver file." << std::endl;
      exit(-1);
    }
  }
  
  return locations;
}

void writeReceiverFile(KDTree const& tree, std::string const& fileName) {
  std::ofstream out(fileName.c_str());
  out << std::scientific << std::setprecision(10);
  
  int failureCounter = 0;
  Point const* points = tree.points();
  for (unsigned p = 0; p < tree.numPoints(); ++p) {
    Point const* point = &points[tree.index(p)];
    if (!std::isnan(point->z)) {
      out << point->x << " " << point->y << " " << point->z << std::endl;
    } else {
      std::cerr << "Warning: Did not find elevation for receiver at (" << point->x << ", " << point->y << ")." << std::endl;
      ++failureCounter;
    }
  }

  if (failureCounter > 0) {
    std::cerr << failureCounter << " points were not written due to missing elevation." << std::endl;
  }
}

void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(-1);
  }
}

Mesh::Mesh(std::string const& fileName){
  int stat;
	size_t maxVertices;
	size_t maxElements;

  stat = nc_open(fileName.c_str(), 0, &ncid);
  check_err(stat,__LINE__,__FILE__);

	int ncDimPart;
	stat = nc_inq_dimid(ncid, "partitions", &ncDimPart);
  check_err(stat,__LINE__,__FILE__);
	stat = nc_inq_dimlen(ncid, ncDimPart, &partitions);
  check_err(stat,__LINE__,__FILE__);

	int ncDimVrtx;
	stat = nc_inq_dimid(ncid, "vertices", &ncDimVrtx);
  check_err(stat,__LINE__,__FILE__);
	stat = nc_inq_dimlen(ncid, ncDimVrtx, &maxVertices);
  check_err(stat,__LINE__,__FILE__);

	int ncDimElem;
	stat = nc_inq_dimid(ncid, "elements", &ncDimElem);
  check_err(stat,__LINE__,__FILE__);
	stat = nc_inq_dimlen(ncid, ncDimElem, &maxElements);
  check_err(stat,__LINE__,__FILE__);

  int ncVarElemSize;
	stat = nc_inq_varid(ncid, "element_size", &ncVarElemSize);
  check_err(stat,__LINE__,__FILE__);

	stat = nc_inq_varid(ncid, "element_vertices", &ncVarElemVertices);
  check_err(stat,__LINE__,__FILE__);

	stat = nc_inq_varid(ncid, "element_boundaries", &ncVarElemBoundaries);
  check_err(stat,__LINE__,__FILE__);

  int ncVarVrtxSize;
	stat = nc_inq_varid(ncid, "vertex_size", &ncVarVrtxSize);
  check_err(stat,__LINE__,__FILE__);

	stat = nc_inq_varid(ncid, "vertex_coordinates", &ncVarVrtxCoords);
  check_err(stat,__LINE__,__FILE__);

	elementSize = new int[partitions];
	stat = nc_get_var_int(ncid, ncVarElemSize, elementSize);
  check_err(stat,__LINE__,__FILE__);

	vertexSize = new int[partitions];
	stat = nc_get_var_int(ncid, ncVarVrtxSize, vertexSize);
  check_err(stat,__LINE__,__FILE__);
  
  vertexCoordinates = new double[3*maxVertices];
  elementBoundaries = new int[4*maxElements];
  elementVertices = new int[4*maxElements];
}

Mesh::~Mesh() {
  int stat = nc_close(ncid);
  check_err(stat,__LINE__,__FILE__);
  
  delete[] elementSize;
  delete[] vertexSize;
  delete[] elementBoundaries;
  delete[] elementVertices;
  delete[] vertexCoordinates;
}

void Mesh::readPartition(int partition) {
  int stat;
  size_t elementsStart[3] = {partition, 0, 0};
  size_t elementsSize[3] = {1, elementSize[partition], 4};
  
  stat = nc_get_vara_int(ncid, ncVarElemVertices, elementsStart, elementsSize, elementVertices);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_get_vara_int(ncid, ncVarElemBoundaries, elementsStart, elementsSize, elementBoundaries);
  check_err(stat,__LINE__,__FILE__);
  
  
  size_t verticesStart[3] = {partition, 0, 0};
  size_t verticesSize[3] = {1, vertexSize[partition], 3};
  
  stat = nc_get_vara_double(ncid, ncVarVrtxCoords, verticesStart, verticesSize, vertexCoordinates);
  check_err(stat,__LINE__,__FILE__);  
}
