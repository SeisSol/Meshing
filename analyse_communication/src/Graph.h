#ifndef GRAPH_H_
#define GRAPH_H_

typedef unsigned long long counter_t;

class Graph {  
private:
  counter_t* m_edges;
  unsigned m_numberOfNodes;

public:
  explicit Graph(unsigned numberOfNodes);
  ~Graph();
  
  void addEdge(unsigned fromNode, unsigned toNode);
  void printDOT(char const* filename);
};

#endif // GRAPH_H_
