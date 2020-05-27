#ifndef _GRAPH_H_
#define _GRAPH_H_

typedef struct {
  unsigned s;
  unsigned t;
} Edge;

typedef struct {
  unsigned n; // Number of nodes
  unsigned e; // Number of edges
  Edge *edges; // List of edges
  unsigned *rank; // Ranking of the nodes according to degeneracy ordering
} EdgeList;

typedef struct {
  unsigned n; // Number of nodes
  unsigned e; // Number of edges
  Edge *edges; // List of edges
  unsigned *cd; // Cumulative degree: cd[0] = 0 and for i >= 1, cd[i] is the sum of degrees of node 0, 1, ..., i - 1
  unsigned *adj; // Truncated list of neighbors
  unsigned core; // Core values of the graph
} Graph;

typedef struct {
  unsigned *n; // n[l]: number of nodes in G_l
  unsigned **d; // d[l]: degrees of ndoes of G_l
  unsigned *adj; // Truncated list of neighbors
  unsigned char *label; // label[i]: label of node i
  unsigned **nodes; // nodes[l]: nodes in G_l
  unsigned core;
} Subgraph;

void FreeGraph(Graph *g);

void FreeSubgraph(Subgraph *sg, unsigned char k);

EdgeList *ReadEdgeList(char *file_name);

// Compute degeneracy ordering.
// The larger the node's core number, the higher (smaller) the node's rank.
void SortByCore(EdgeList *el);

// Relabel each node to be its rank in degeneracy ordering.
// Each edge is oriented from the endpoint with lower (larger) rank to that with higher (smaller) rank.
void Relabel(EdgeList *el);

// Building the directed graph.
Graph *MakeGraph(EdgeList *el);

// Count the number of edges in a node-reduced subgraph
unsigned CountEdges(Graph *g, unsigned num_nodes, unsigned nodes[]);

#endif
