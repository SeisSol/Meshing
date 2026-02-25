#include <iostream>
#include <string>

#include <mpi.h>

#include "mesh.h"
#include "utils/args.h"
#include "utils/logger.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  utils::Args args;
  args.addOption("mesh", 'm', "mesh file (h5)");

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  const auto meshFile = args.getArgument<std::string>("mesh");

  if (rank == 0) {
    logInfo() << "Read" << meshFile << ".";
  }
  const Mesh mesh(meshFile);
  const bool correctOnRank = mesh.checkNeighbors();
  bool correctOnAllRanks{true};

  MPI_Reduce(&correctOnRank, &correctOnAllRanks, 1, MPI_CXX_BOOL, MPI_LAND, 0, MPI_COMM_WORLD);
  MPI_Finalize();

  if (rank == 0) {
    if (correctOnAllRanks) {
      logInfo() << "Mesh is correct";
      return 0;
    } else {
      logInfo() << "Found at least one broken element, see above for more details.";
      return 1;
    }
  }
  return 0;
}
