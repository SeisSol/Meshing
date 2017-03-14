#include <iostream>
#include <utils/args.h>
#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/random.hpp>

#include "common/PartitionReader.h"

double calculateTimestep(double cfl, double insphereRadius, double pVelocity, unsigned order)
{
  return cfl * 2.0 * insphereRadius / (pVelocity) / (2*(order-1)+1);
}

unsigned getCluster(double timestep, double globalMinTimestep, unsigned rate)
{
  if (rate == 1) {
    return 0;
  }

  double upper;
  upper = rate * globalMinTimestep;
  
  unsigned cluster = 0;
  while (upper <= timestep) {
    upper *= rate;
    ++cluster;
  }
  return cluster;
}

int main(int argc, char** argv)
{
  utils::Args args;
  args.addOption("mesh", 'm', "Netcdf mesh file");
  args.addOption("average-flops-per-cell", 'f', "Hardware flops per cell");
  args.addOption("node-performance", 'n', "Performance of a node in GFLOPS");
  args.addOption("p-wave-velocity", 'p', "P-wave velocity (m/s)", utils::Args::Required, false);
  args.addOption("CFL", 'c', "CFL number", utils::Args::Required, false);
  args.addOption("order", 'o', "Convergence order", utils::Args::Required, false);
  args.addOption("rate", 'r', "Clusterd lts rate", utils::Args::Required, false);
  args.addOption("final-time", 't', "Final simulation time", utils::Args::Required, false);

	if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }
  
  std::string meshFile = args.getArgument<std::string>("mesh");  
  double flopsPerCell = args.getArgument<double>("average-flops-per-cell");
  double nodePerformance = args.getArgument<double>("node-performance");
  double pVelocity = args.getArgument<double>("p-wave-velocity", 6000.0);
  double cfl = args.getArgument<double>("CFL", 0.5);
  unsigned order = args.getArgument<unsigned>("order", 6);
  unsigned rate = args.getArgument<unsigned>("rate", 2);
  double finalTime = args.getArgument<unsigned>("final-time", 1.0);
  
  PartitionReader reader(meshFile);
  std::vector<double> timesteps;

  unsigned numberOfElements = 0;
  double minInsphereRadius = std::numeric_limits<double>::max();
  double minTimestep = std::numeric_limits<double>::max();
  double maxTimestep = std::numeric_limits<double>::min();
  for (unsigned p = 0; p < reader.partitions; ++p) {    
    reader.readPartition(p);
    numberOfElements += reader.elementSize[p];

    for (unsigned elem = 0; elem < reader.elementSize[p]; ++elem) {
      glm::dvec3 x[4];
      for (unsigned vtx = 0; vtx < 4; ++vtx) {
        for (unsigned d = 0; d < 3; ++d) {
          x[vtx][d] = reader.vertexCoordinates[3 * reader.elementVertices[4*elem + vtx] + d];
        }
      }

      double alpha = determinant(glm::dmat4(glm::dvec4(x[0], 1.0), glm::dvec4(x[1], 1.0), glm::dvec4(x[2], 1.0), glm::dvec4(x[3], 1.0)));
      double Nabc = length(cross(x[1]-x[0], x[2]-x[0]));
      double Nabd = length(cross(x[1]-x[0], x[3]-x[0]));
      double Nacd = length(cross(x[2]-x[0], x[3]-x[0]));
      double Nbcd = length(cross(x[2]-x[1], x[3]-x[1]));
      double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);
      minInsphereRadius = std::min(minInsphereRadius, insphere);
      
      double timestep = calculateTimestep(cfl, insphere, pVelocity, order);
      minTimestep = std::min(minTimestep, timestep);
      maxTimestep = std::max(maxTimestep, timestep);
      
      timesteps.push_back(timestep);
    }
  }
  
  unsigned numClusters = getCluster(maxTimestep, minTimestep, rate);
  std::vector<unsigned> clusterHistogram(numClusters+1, 0);
  for (std::vector<double>::const_iterator it = timesteps.begin(); it != timesteps.end(); ++it) {
    unsigned cluster = getCluster(*it, minTimestep, rate);
    ++clusterHistogram[cluster];
  }
  
  unsigned numberOfUpdates = 0;
  double clusterTimestep = minTimestep;
  for (std::vector<unsigned>::const_iterator it = clusterHistogram.begin(); it != clusterHistogram.end(); ++it) {
    numberOfUpdates += (*it) * std::ceil(finalTime / clusterTimestep);
    clusterTimestep *= rate;
  }
  
  std::cout << "Minimum insphere radius: " << minInsphereRadius << std::endl;
  std::cout << "Elements: " << numberOfElements << std::endl;
  std::cout << "Minimum timestep: " << minTimestep << std::endl;
  std::cout << "Maximum timestep: " << maxTimestep << std::endl;
  std::cout << "Elements in time clusters: ";
  for (std::vector<unsigned>::const_iterator it = clusterHistogram.begin(); it != clusterHistogram.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
  std::cout << "Number of cell updates (@ " << finalTime << "s): " << numberOfUpdates << std::endl;
  std::cout << "Estimated time on one node: " << numberOfUpdates * flopsPerCell / (nodePerformance * 1.0e9) << std::endl;
  
  return 0;
}
