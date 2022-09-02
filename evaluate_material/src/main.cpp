#include <iostream>
#include <string>

#include <mpi.h>
#include <utils/args.h>
#include <utils/logger.h>

#include "mesh.h"
#include "parameterDB.h"
#include "writer.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  utils::Args args;
  args.addOption("mesh", 'm', "mesh file (h5)");
  args.addOption("easi", 'e', "easi file (yaml)");
  args.addOption("output", 'o', "output file prefix");

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  const auto meshFile = args.getArgument<std::string>("mesh");
  const auto easiFile = args.getArgument<std::string>("easi");
  const auto outputFilePrefix = args.getArgument<std::string>("output");

  logInfo() << "Read" << meshFile << ".";
  const Mesh mesh(meshFile);
  const auto barycenters = mesh.getElementBarycenters();
  logInfo() << "Done";
  logInfo() << "Evaluate material from" << easiFile << "on the mesh.";
  ParameterDB parameterDB(easiFile);
  auto query = parameterDB.generateQuery(barycenters, mesh.getElementGroups());
  const auto [parameters, material] = parameterDB.evaluate(query);
  logInfo() << "Done";
  logInfo() << "Write result to " << outputFilePrefix << ".";
  Writer writer(outputFilePrefix, mesh);
  writer.write(parameters, material);
  logInfo() << "Done";

  MPI_Finalize();

  return 0;
}
