#ifndef PARTITIONREADER_H_
#define PARTITIONREADER_H_

#include <string>

class PartitionReader {
public:
  explicit PartitionReader(std::string const& fileName);
  ~PartitionReader();
  
  void readPartition(int partition);

	size_t partitions;
  int* elementSize;
  int* vertexSize;
  int* elementBoundaries;
  int* elementVertices;
  int* elementNeighborRanks;
  double* vertexCoordinates;
  
private:
  int ncid;
	int ncVarElemVertices;
	int ncVarElemBoundaries;
	int ncVarElemNeighborRanks;
	int ncVarVrtxCoords;
};

#endif // PARTITIONREADER_H_
