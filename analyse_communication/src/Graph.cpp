#include "Graph.h"

#include <cstdio>
#include <cstring>
#include <algorithm>

Graph::Graph(unsigned numberOfNodes)
  : m_numberOfNodes(numberOfNodes)
{
  m_edges = new counter_t[numberOfNodes * numberOfNodes];
  memset(m_edges, 0, numberOfNodes * numberOfNodes * sizeof(counter_t));
}

Graph::~Graph()
{
  delete[] m_edges;
}

void Graph::addEdge(unsigned fromNode, unsigned toNode)
{
  ++m_edges[fromNode * m_numberOfNodes + toNode];
}

void Graph::printDOT(char const* filename)
{
  FILE* fp = fopen(filename, "w");
  fprintf(fp, "graph communication {\n");
  fprintf(fp, "overlap=scale;\n");
  fprintf(fp, "splines=true;\n");
  counter_t maxEdgeCut = 0;
  for (unsigned from = 0; from < m_numberOfNodes; ++from) {
    for (unsigned to = from+1; to < m_numberOfNodes; ++to) {
      maxEdgeCut = std::max(maxEdgeCut, m_edges[from * m_numberOfNodes + to] + m_edges[to * m_numberOfNodes + from]);
    }
  }
  for (unsigned from = 0; from < m_numberOfNodes; ++from) {
    counter_t edgeCut = 0;
    for (unsigned to = 0; to < m_numberOfNodes; ++to) {
      edgeCut += m_edges[from * m_numberOfNodes + to];
    }
    fprintf(fp, "%u [label=\"P%u\n%u\"];\n", from, from, edgeCut);
    
    for (unsigned to = from+1; to < m_numberOfNodes; ++to) {
      counter_t edgeCut = m_edges[from * m_numberOfNodes + to] + m_edges[to * m_numberOfNodes + from];
      if (edgeCut > 0) {
        double edgeWeight = static_cast<double>(edgeCut) / maxEdgeCut;
        fprintf(fp, "%u -- %u [weight=%lf,penwidth=%lf];\n", from, to, edgeWeight, 7.0*edgeWeight);
      }
    }
  }
 
  fprintf(fp, "}\n");
  
  fclose(fp);
}
