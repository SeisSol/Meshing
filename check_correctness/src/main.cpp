#include <iostream>
#include <string>

#include <mpi.h>

#include "mesh.h"
#include "utils/args.h"
#include "utils/logger.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  utils::Args args;
  args.addOption("mesh", 'm', "mesh file (h5)");

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  const auto meshFile = args.getArgument<std::string>("mesh");

  logInfo() << "Read" << meshFile << ".";
  const Mesh mesh(meshFile);
  mesh.checkNeighbors();

  MPI_Finalize();

  return 0;
}
