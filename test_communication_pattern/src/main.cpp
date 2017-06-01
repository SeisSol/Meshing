#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <utils/args.h>
#include <utils/logger.h>

#include "Comm.h"

counter_t* readCommMatrix(std::string const& matrixFile)
{
  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  
  counter_t* edges = NULL;
  FILE* file = fopen(matrixFile.c_str(), "rb");
  fseek(file, 0, SEEK_END);
  long int size = ftell(file);
  
  int N = static_cast<int>(sqrt(size / sizeof(counter_t)));

  if (N*N*sizeof(counter_t) != size) {
    logError() << "The communication matrix is not a square matrix.";
    return NULL;
  }  
  if (N != numRanks) {
    logError() << "Communication matrix has" << N << "ranks whereas this test was started with" << numRanks << "ranks.";
    return NULL;
  }

  fseek(file, rank * N * sizeof(counter_t), SEEK_SET);
  edges = new counter_t[N];
  fread(edges, sizeof(counter_t), N, file);
  fclose(file);
  
  return edges;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
    
  utils::Args args;
  args.addOption("order", 'o', "Convergence order.", utils::Args::Required, false);
  args.addOption("num-quantities", 'q', "Number of quantities.", utils::Args::Required, false);
  args.addAdditionalOption("matrix", "Edge-cut matrix.");

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  int order = args.getArgument<int>("order", 6);
  int nq = args.getArgument<int>("num-quantities", 9);
  std::string matrixFile = args.getAdditionalArgument<std::string>("matrix");
  
  counter_t* edges = readCommMatrix(matrixFile);
  
  size_t messageSize = sizeof(double) * nq * order * (order+1) * (order+2) / 6;
  testCommunication(edges, messageSize);
  
  delete[] edges;  
  
  MPI_Finalize();
  
  return 0;
}
