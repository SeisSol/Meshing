#include <iostream>
#include <utils/args.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "SeismicVelocity.h"
#include "Writer.h"
#include "common/PartitionReader.h"

double calculateTimestep(double cfl, double insphereRadius, double pVelocity, unsigned order)
{
  return cfl * 2.0 * insphereRadius / (pVelocity) / (2*(order-1)+1);
}

int main(int argc, char** argv)
{
  utils::Args args;
  args.addOption("mesh", 'm', "Netcdf mesh file");
  args.addOption("velocity-model", 'v', "Velocity model");
  args.addOption("CFL", 'c', "CFL number", utils::Args::Required, false);
  args.addOption("order", 'O', "Convergence order");
  args.addOption("output", 'o', "Output file name");

	if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }
  
  std::string meshFile = args.getArgument<std::string>("mesh");
  std::string velocityModel = args.getArgument<std::string>("velocity-model");
  double cfl = args.getArgument<double>("CFL", 0.5);
  unsigned order = args.getArgument<unsigned>("order");
  std::string outputFile = args.getArgument<std::string>("output");
  
  double (*pWaveVelocityFun)(int,double,double,double);
  if (velocityModel == "landers61") {
    pWaveVelocityFun = &landers61;
  } else if (velocityModel == "sumatra1223_low") {
    pWaveVelocityFun = &sumatra1223_low;
  } else if (velocityModel == "sumatra1223_high") {
    pWaveVelocityFun = &sumatra1223_high;
  } else if (velocityModel == "sumatra1224") {
    pWaveVelocityFun = &sumatra1224;
  } else {
    std::cerr << "Error: Unknown velocity model." << std::endl;
    exit(-1);
  }
  
  PartitionReader reader(meshFile);
  std::vector<double> timesteps;

  unsigned numberOfElements = 0;
  for (unsigned p = 0; p < reader.partitions; ++p) {
    std::cout << "Reading partition " << p+1 << "..." << std::endl; 
    reader.readPartition(p);
    numberOfElements += reader.elementSize[p];

    for (unsigned elem = 0; elem < reader.elementSize[p]; ++elem) {
      glm::dvec3 barycentre(0.,0.,0.);
      glm::dvec3 x[4];
      for (unsigned vtx = 0; vtx < 4; ++vtx) {
        for (unsigned d = 0; d < 3; ++d) {
          x[vtx][d] = reader.vertexCoordinates[3 * reader.elementVertices[4*elem + vtx] + d];
        }
        barycentre += 0.25 * x[vtx];
      }

      double alpha = determinant(glm::dmat4(glm::dvec4(x[0], 1.0), glm::dvec4(x[1], 1.0), glm::dvec4(x[2], 1.0), glm::dvec4(x[3], 1.0)));
      double Nabc = length(cross(x[1]-x[0], x[2]-x[0]));
      double Nabd = length(cross(x[1]-x[0], x[3]-x[0]));
      double Nacd = length(cross(x[2]-x[0], x[3]-x[0]));
      double Nbcd = length(cross(x[2]-x[1], x[3]-x[1]));
      double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);
      
      double timestep = calculateTimestep(cfl, insphere, pWaveVelocityFun(reader.elementGroup[elem], barycentre.x, barycentre.y, barycentre.z), order);
      
      timesteps.push_back(timestep);
    }
  }
  
  writeTimesteps(timesteps, outputFile);
  
  return 0;
}
