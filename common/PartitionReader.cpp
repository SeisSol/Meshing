/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
#include "PartitionReader.h"

#include <cstdlib>
#include <cstdio>
#include <netcdf.h>

void check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(-1);
  }
}

PartitionReader::PartitionReader(std::string const& fileName){
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

	stat = nc_inq_varid(ncid, "element_neighbor_ranks", &ncVarElemNeighborRanks);
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
  elementNeighborRanks = new int[4*maxElements];
}

PartitionReader::~PartitionReader() {
  int stat = nc_close(ncid);
  check_err(stat,__LINE__,__FILE__);
  
  delete[] elementSize;
  delete[] vertexSize;
  delete[] elementBoundaries;
  delete[] elementVertices;
  delete[] elementNeighborRanks;
  delete[] vertexCoordinates;
}

void PartitionReader::readPartition(int partition) {
  int stat;
  size_t elementsStart[3] = {partition, 0, 0};
  size_t elementsSize[3] = {1, elementSize[partition], 4};
  
  stat = nc_get_vara_int(ncid, ncVarElemVertices, elementsStart, elementsSize, elementVertices);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_get_vara_int(ncid, ncVarElemBoundaries, elementsStart, elementsSize, elementBoundaries);
  check_err(stat,__LINE__,__FILE__);
  
  stat = nc_get_vara_int(ncid, ncVarElemNeighborRanks, elementsStart, elementsSize, elementNeighborRanks);
  check_err(stat,__LINE__,__FILE__);
  
  
  size_t verticesStart[3] = {partition, 0, 0};
  size_t verticesSize[3] = {1, vertexSize[partition], 3};
  
  stat = nc_get_vara_double(ncid, ncVarVrtxCoords, verticesStart, verticesSize, vertexCoordinates);
  check_err(stat,__LINE__,__FILE__);  
}
