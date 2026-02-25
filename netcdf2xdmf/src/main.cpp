/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <valarray>

#include <netcdf.h>

#include "utils/args.h"
#include "utils/logger.h"
#include "utils/path.h"
#include "utils/stringutils.h"

#include "easi/YAMLParser.h"

int main(int argc, char* argv[])
{
	utils::Args args;
	args.addOption("boundary", 'b', "Convert only boundary cells", utils::Args::No, false);
	args.addOption("materialFile", 'm', "Add material parameters to xdmf", utils::Args::Required, false);
	args.addAdditionalOption("input", "The netCDF mesh file");
	args.addAdditionalOption("output", "The generated XDMF file", false);

	switch (args.parse(argc, argv)) {
	case utils::Args::Error:
		return 1;
	case utils::Args::Help:
		return 2;
	}

	bool convertBoundaries = args.isSet("boundary");
	std::string material = args.getArgument<std::string>("materialFile", "");

	// Get infile/outfile
	std::string infile = args.getAdditionalArgument<std::string>("input");
	std::string outfilePrefix;
	if (args.isSetAdditional("output")) {
		outfilePrefix = args.getAdditionalArgument<std::string>("output");
	} else {
		outfilePrefix = infile;
		size_t pos = outfilePrefix.find_last_of('.');
		outfilePrefix.resize(pos);
	}

	// Remove file ending
	if (utils::StringUtils::endsWith(outfilePrefix, ".xdmf")) {
		outfilePrefix.resize(outfilePrefix.size()-5);
	} else {
		if (convertBoundaries) {
			if (!utils::StringUtils::endsWith(outfilePrefix, "_bnd"))
				outfilePrefix += "_bnd";
		}
	}

	std::string outfilePrefixShort = utils::Path(outfilePrefix).basename();

	// Open the xml file
	std::ofstream xdmfFile((outfilePrefix+".xdmf").c_str());

	xdmfFile << "<?xml version=\"1.0\" ?>" << std::endl
		<< "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
		<< "<Xdmf Version=\"2.0\">" << std::endl
		<< " <Domain>" << std::endl;
	xdmfFile << "  <Grid Name=\"netcdf\" GridType=\"Uniform\">" << std::endl;

//
// 			m_xdmfFile << "  <DataItem NumberType=\"UInt\" Precision=\"4\" Format=\""
// 						<< backends::Backend::format() << "\" Dimensions=\"" << totalSize[0] << "\">"
// 					<< backends::Backend::dataItemLocation(backendPrefix.c_str(), "partition")
// 					<< "</DataItem>" << std::endl;

	// Open the netCDF file
	int ncFile;
	nc_open(infile.c_str(), 0, &ncFile);

	int ncDimPart;
	nc_inq_dimid(ncFile, "partitions", &ncDimPart);
	size_t partitions;
	nc_inq_dimlen(ncFile, ncDimPart, &partitions);

	int ncDimVrtx;
	nc_inq_dimid(ncFile, "vertices", &ncDimVrtx);
	size_t maxVertices;
	nc_inq_dimlen(ncFile, ncDimVrtx, &maxVertices);

	int ncDimElem;
	nc_inq_dimid(ncFile, "elements", &ncDimElem);
	size_t maxElements;
	nc_inq_dimlen(ncFile, ncDimElem, &maxElements);

	// Create netcdf variables
	int ncVarElemSize;
	nc_inq_varid(ncFile, "element_size", &ncVarElemSize);

	int ncVarElemVertices;
	nc_inq_varid(ncFile, "element_vertices", &ncVarElemVertices);

	int ncVarElemBoundaries;
	nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries);

	int ncVarElemGroup;
	bool hasGroup = false;
	int ncResult = nc_inq_varid(ncFile, "element_group", &ncVarElemGroup);
	if (ncResult != NC_ENOTVAR)
		hasGroup = true;

	int ncVarVrtxSize;
	nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize);

	int ncVarVrtxCoords;
	nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords);

	unsigned int* vertexSize = new unsigned int[partitions];
	nc_get_var_uint(ncFile, ncVarVrtxSize, vertexSize);
	unsigned int* elementSize = new unsigned int[partitions];
	nc_get_var_uint(ncFile, ncVarElemSize, elementSize);

	// Write cells
	std::string topoName;
	unsigned int topoSize;

	unsigned int nvertices = std::valarray<unsigned int>(vertexSize, partitions).sum();
	unsigned int nelements = 0;
	if (convertBoundaries) {
		topoName = "Triangle";
		topoSize = 3;

		unsigned int* boundaries = new unsigned int[maxElements*4];

		size_t start[3] = {0, 0, 0};
		size_t size[3] = {1, 0, 4};
		for (size_t i = 0; i < partitions; i++) {
			start[0] = i;
			size[1] = elementSize[i];
			nc_get_vara_uint(ncFile, ncVarElemBoundaries, start, size, boundaries);

			for (unsigned int j = 0; j < elementSize[i]; j++) {
				for (unsigned int k = 0; k < 4; k++) {
					if (boundaries[j*4+k] != 0)
						nelements++;
				}
			}
		}

		delete [] boundaries;
	} else {
		topoName = "Tetrahedron";
		topoSize = 4;

		nelements = std::valarray<unsigned int>(elementSize, partitions).sum();
	}

	xdmfFile << "   <Topology TopologyType=\"" << topoName << "\" NumberOfElements=\"" << nelements << "\">" << std::endl
		// This should be UInt but for some reason this does not work with binary data
		<< "    <DataItem NumberType=\"Int\" Precision=\"4\" Format=\"Binary\" Dimensions=\"" << nelements << " " << topoSize << "\">"
		<< outfilePrefixShort << "_connect.bin</DataItem>" << std::endl
		<< "   </Topology>" << std::endl;

	int fd = open((outfilePrefix+"_connect.bin").c_str(), O_CREAT | O_WRONLY | O_TRUNC,
		S_IRUSR |S_IWUSR |S_IRGRP |S_IWGRP | S_IROTH | S_IWOTH);

	int* elements = new int[maxElements*4];
	unsigned int vrtxIdStart = 0;
	if (convertBoundaries) {
		int* connect = new int[nelements*3];
		unsigned int nextElement = 0;

		unsigned int* boundaries = new unsigned int[maxElements*4];

		for (size_t i = 0; i < partitions; i++) {
			size_t start[3] = {i, 0, 0};
			size_t size[3] = {1, elementSize[i], 4};
			nc_get_vara_int(ncFile, ncVarElemVertices, start, size, elements);
			nc_get_vara_uint(ncFile, ncVarElemBoundaries, start, size, boundaries);

			for (unsigned int j = 0; j < elementSize[i]; j++) {
				if (boundaries[j*4] != 0) {
					connect[nextElement*3] = vrtxIdStart+elements[j*4];
					connect[nextElement*3+1] = vrtxIdStart+elements[j*4+1];
					connect[nextElement*3+2] = vrtxIdStart+elements[j*4+2];
					nextElement++;
				}
				if (boundaries[j*4+1] != 0) {
					connect[nextElement*3] = vrtxIdStart+elements[j*4];
					connect[nextElement*3+1] = vrtxIdStart+elements[j*4+1];
					connect[nextElement*3+2] = vrtxIdStart+elements[j*4+3];
					nextElement++;
				}
				if (boundaries[j*4+2] != 0) {
					connect[nextElement*3] = vrtxIdStart+elements[j*4];
					connect[nextElement*3+1] = vrtxIdStart+elements[j*4+2];
					connect[nextElement*3+2] = vrtxIdStart+elements[j*4+3];
					nextElement++;
				}
				if (boundaries[j*4+3] != 0) {
					connect[nextElement*3] = vrtxIdStart+elements[j*4+1];
					connect[nextElement*3+1] = vrtxIdStart+elements[j*4+2];
					connect[nextElement*3+2] = vrtxIdStart+elements[j*4+3];
					nextElement++;
				}

			}

			vrtxIdStart += vertexSize[i];
		}

		delete [] boundaries;

		write(fd, connect, nelements*3*sizeof(int));
	} else {
		for (size_t i = 0; i < partitions; i++) {
			size_t start[3] = {i, 0, 0};
			size_t size[3] = {1, elementSize[i], 4};
			nc_get_vara_int(ncFile, ncVarElemVertices, start, size, elements);

			for (unsigned int j = 0; j < elementSize[i]*4; j++)
				elements[j] = vrtxIdStart + elements[j];

			vrtxIdStart += vertexSize[i];

			write(fd, elements, elementSize[i]*4*sizeof(int));
		}
	}
	delete [] elements;

	close(fd);

	// Write vertices
	xdmfFile << "   <Geometry name=\"geo\" GeometryType=\"XYZ\" NumberOfElements=\"" << nvertices << "\">" << std::endl
		<< "    <DataItem NumberType=\"Float\" Precision=\"4\" Format=\"Binary\" Dimensions=\"" << nvertices << " 3\">"
		<< outfilePrefixShort << "_geometry.bin</DataItem>" << std::endl
		<< "   </Geometry>" << std::endl;

	fd = open((outfilePrefix+"_geometry.bin").c_str(), O_CREAT | O_WRONLY | O_TRUNC,
		S_IRUSR |S_IWUSR |S_IRGRP |S_IWGRP | S_IROTH | S_IWOTH);

	float* vertices = new float[maxVertices*3];
	for (size_t i = 0; i < partitions; i++) {
		size_t start[3] = {i, 0, 0};
		size_t size[3] = {1, vertexSize[i], 3};
		nc_get_vara_float(ncFile, ncVarVrtxCoords, start, size, vertices);

		write(fd, vertices, vertexSize[i]*3*sizeof(float));
	}

	delete [] vertices;

	close(fd);

	// Write Partition/Boundary type
	std::string  attributeName;
	if (convertBoundaries)
		attributeName = "boundary";
	else
		attributeName = "partition";

	xdmfFile << "   <Attribute Name=\"" << attributeName << "\" Center=\"Cell\">" << std::endl
		<< "    <DataItem NumberType=\"Int\" Precision=\"4\" Format=\"Binary\" Dimensions=\"" << nelements << "\">"
		<< outfilePrefixShort << "_" << attributeName << ".bin</DataItem>" << std::endl
		<< "   </Attribute>" << std::endl;

	fd = open((outfilePrefix+"_"+attributeName+".bin").c_str(), O_CREAT | O_WRONLY | O_TRUNC,
		S_IRUSR |S_IWUSR |S_IRGRP |S_IWGRP | S_IROTH | S_IWOTH);

	if (convertBoundaries) {
		unsigned int* boundaries = new unsigned int[maxElements*4];

		int* boundariesOut = new int[nelements];
		unsigned int nextElement = 0;

		for (size_t i = 0; i < partitions; i++) {
			size_t start[3] = {i, 0, 0};
			size_t size[3] = {1, elementSize[i], 4};
			nc_get_vara_uint(ncFile, ncVarElemBoundaries, start, size, boundaries);

			for (unsigned int j = 0; j < elementSize[i]; j++) {
				for (unsigned int k = 0; k < 4; k++) {
					if (boundaries[j*4+k] != 0)
						boundariesOut[nextElement++] = boundaries[j*4+k];
				}
			}
		}

		write(fd, boundariesOut, nelements*sizeof(int));

		delete [] boundaries;
		delete [] boundariesOut;
	} else {
		int* partition = new int[maxElements];

		for (size_t i = 0; i < partitions; i++) {
			std::fill_n(partition, elementSize[i], i);

			write(fd, partition, elementSize[i]*sizeof(int));
		}

		// Write groups if available
		if (hasGroup) {
			close(fd);

			xdmfFile << "   <Attribute Name=\"group\" Center=\"Cell\">" << std::endl
				<< "    <DataItem NumberType=\"Int\" Precision=\"4\" Format=\"Binary\" Dimensions=\"" << nelements << "\">"
				<< outfilePrefixShort << "_group.bin</DataItem>" << std::endl
				<< "   </Attribute>" << std::endl;

			fd = open((outfilePrefix+"_group.bin").c_str(), O_CREAT | O_WRONLY | O_TRUNC,
				S_IRUSR |S_IWUSR |S_IRGRP |S_IWGRP | S_IROTH | S_IWOTH);

			int* groups = new int[maxElements];
			size_t start[3] = {0, 0};
			size_t size[3] = {1, 0};
			for (size_t i = 0; i < partitions; i++) {
				start[0] = i;
				size[1] = elementSize[i];
				nc_get_vara_int(ncFile, ncVarElemGroup, start, size, groups);

				write(fd, groups, elementSize[i]*sizeof(int));
			}

			delete [] groups;
		}
    
    if (material.size() > 0) {     
      easi::YAMLParser parser(3);
      easi::Component* model = parser.parse(material);
      
      std::vector<std::string> paramNames{"lambda", "mu", "rho"};    

      int fds[3];
      for (unsigned d = 0; d < 3; ++d) {
        xdmfFile << "   <Attribute Name=\"" << paramNames[d] << "\" Center=\"Cell\">" << std::endl
          << "    <DataItem NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Dimensions=\"" << nelements << "\">"
          << outfilePrefixShort << "_" << paramNames[d] << ".bin</DataItem>" << std::endl
          << "   </Attribute>" << std::endl;

        fds[d] = open((outfilePrefix+"_" + paramNames[d] + ".bin").c_str(), O_CREAT | O_WRONLY | O_TRUNC,
          S_IRUSR |S_IWUSR |S_IRGRP |S_IWGRP | S_IROTH | S_IWOTH);
      }
      
      easi::ArraysAdapter adapter;
      double* parameters[3];
      for (unsigned d = 0; d < 3; ++d) {
        parameters[d] = new double[maxElements];
        adapter.addBindingPoint(paramNames[d], parameters[d]);
      }
      
      int* elements = new int[maxElements*4];
      float* vertices = new float[maxVertices*3];
      int* groups;
      if (hasGroup) {
        groups = new int[maxElements];
      }
      for (size_t i = 0; i < partitions; ++i) {
        size_t startV[3] = {i, 0, 0}; size_t sizeV[3] = {1, elementSize[i], 4};
        nc_get_vara_int(ncFile, ncVarElemVertices, startV, sizeV, elements);
        size_t startE[3] = {i, 0, 0}; size_t sizeE[3] = {1, vertexSize[i], 3};
        nc_get_vara_float(ncFile, ncVarVrtxCoords, startE, sizeE, vertices);
        if (hasGroup) {
          size_t startG[2] = {i, 0}; size_t sizeG[2] = {1, elementSize[i]};
          nc_get_vara_int(ncFile, ncVarElemGroup, startG, sizeG, groups);
        }

        easi::Query query(elementSize[i], 3);
        if (hasGroup) {
          for (unsigned int j = 0; j < elementSize[i]; ++j) {
            query.group(j) = groups[j];
          }          
        }
        for (unsigned int j = 0; j < elementSize[i]; ++j) {
          for (unsigned d = 0; d < 3; ++d) {
            query.x(j,d) = 0.25f * vertices[ 3*elements[4*j] + d];
          }
          for (unsigned v = 1; v < 4; ++v) {
            for (unsigned d = 0; d < 3; ++d) {
              query.x(j,d) += 0.25f * vertices[ 3*elements[4*j+v] + d];
            }
          }
        }
        
        model->evaluate(query, adapter);
        
        for (unsigned d = 0; d < 3; ++d) {
          write(fds[d], parameters[d], elementSize[i]*sizeof(double));
        }
      }
      delete [] elements;
      delete [] vertices;
      if (hasGroup) {
        delete [] groups;
      }
      for (unsigned d = 0; d < 3; ++d) {
        delete[] parameters[d];
        close(fds[d]);
      }
      
      delete model;
    }
	}

	close(fd);

	xdmfFile << "  </Grid>" << std::endl
		<< " </Domain>" << std::endl
		<< "</Xdmf>" << std::endl;

	logInfo() << "Finished successfully";
	return 0;
}
