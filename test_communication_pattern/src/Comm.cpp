#include <mpi.h>
#include "Comm.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <sys/time.h>
#include <algorithm>
#include <limits>

inline double usec(struct timeval start, struct timeval end) {
  return ((double)(((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec))));
}

void testCommunication(counter_t* edgeCut, size_t messageSize, int numTimesteps)
{
  struct timeval start_time, end_time, start_time_total, end_time_total;
  double time_min = std::numeric_limits<double>::max(), time_max = std::numeric_limits<double>::min(), time_avg;
  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  
  unsigned N = 0;
  counter_t E = 0;
  for (unsigned r = 0; r < numRanks; ++r) {
    N += (edgeCut[r] > 0) ? 1 : 0;
    E += edgeCut[r];
  }
  
  MPI_Request* requests = new MPI_Request[2*N];
  unsigned char* copy = new unsigned char[E*messageSize];
  unsigned char* ghost = new unsigned char[E*messageSize];
  for (unsigned i = 0; i < E*messageSize; ++i) {
    copy[i] = static_cast<unsigned char>(lrand48() % 256);
    ghost[i] = static_cast<unsigned char>(lrand48() % 256);
  }
  
  if (rank == 0) {
    printf("Rank       tavg       tmin       tmax         bw\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
 
  gettimeofday(&start_time_total, NULL);
  for (unsigned t = 0; t < numTimesteps; ++t) {
    unsigned e = 0;
    unsigned n = 0;
    gettimeofday(&start_time, NULL);
    for (unsigned r = 0; r < numRanks; ++r) {
      if (edgeCut[r] > 0) {
        MPI_Isend(copy + e*messageSize, edgeCut[r]*messageSize, MPI_BYTE, r, 0, MPI_COMM_WORLD, requests + 2*n);
        MPI_Irecv(ghost + e*messageSize, edgeCut[r]*messageSize, MPI_BYTE, r, 0, MPI_COMM_WORLD, requests + 2*n+1);
        e += edgeCut[r];
        ++n;
      }
    }
    MPI_Waitall(2*N, requests, MPI_STATUSES_IGNORE);
    gettimeofday(&end_time, NULL);
    double time = usec(start_time, end_time);
    time_min = std::min(time_min, time);
    time_max = std::max(time_max, time);
  }
  gettimeofday(&end_time_total, NULL);
  time_avg = usec(start_time_total, end_time_total) / numTimesteps;
  
  printf("%4d %10.2lf %10.2lf %10.2lf %10.2lf\n", rank, time_avg, time_min, time_max, E*messageSize / 1.048576 / time_avg);
  
  delete[] requests;    
}
