#include <stdio.h>
#include <stdlib.h>
#include "BinaryHeap.h"
#include "Graph.h"

const unsigned MAX_EDGES = 100000000; // Maximum number of edges for memory allocation; will increase if needed

void FreeGraph(Graph *g) {
  free(g->edges);
  free(g->cd);
  free(g->adj);
  free(g);
}

void FreeSubgraph(Subgraph *sg, unsigned char k) {
  unsigned char i;
  free(sg->n);
  for (i = 1; i < k; ++i){
    free(sg->d[i]);
    free(sg->nodes[i]);
  }
  free(sg->d);
  free(sg->nodes);
  free(sg->label);
  free(sg->adj);
  free(sg);
}

inline unsigned max3(unsigned a, unsigned b, unsigned c) {
  a = (a > b) ? a : b;
  return (a > c) ? a : c;
}

EdgeList *ReadEdgeList(char *file_name) {
  unsigned e1 = MAX_EDGES;
  EdgeList *el = (EdgeList *)malloc(sizeof(EdgeList));
  FILE *file = fopen(file_name, "r");
  el->n = 0;
  el->e = 0;
  el->edges = (Edge *)malloc(e1 * sizeof(Edge));
  while (fscanf(file, "%u %u", &(el->edges[el->e].s), &(el->edges[el->e].t)) == 2) { // Read one edge
    el->n = max3(el->n, el->edges[el->e].s, el->edges[el->e].t);
    ++el->e;
    if (el->e == e1) {
      e1 += MAX_EDGES;
      el->edges = (Edge *)realloc(el->edges, e1 * sizeof(Edge));
    }
  }
  fclose(file);
  ++el->n; // So the nodes are numbered from 0 to el->n - 1
  el->edges = (Edge *)realloc(el->edges, el->e * sizeof(Edge));
  return el;
}

// Compute degeneracy ordering.
// The larger the node's core number, the higher (smaller) the node's rank.
void SortByCore(EdgeList *el) {
  unsigned n = el->n, e = el->e;
  unsigned *d = (unsigned *)calloc(n, sizeof(unsigned));
  unsigned *cd = (unsigned *)malloc((n + 1) * sizeof(unsigned));
  unsigned *adj = (unsigned *)malloc(2 * e * sizeof(unsigned));
  for (unsigned i = 0; i < e; ++i) {
    ++d[el->edges[i].s];
    ++d[el->edges[i].t];
  }
  cd[0] = 0;
  for (unsigned i = 1; i < n + 1; ++i) {
    cd[i] = cd[i - 1] + d[i - 1];
    d[i - 1] = 0;
  }
  for (unsigned i = 0; i < e; ++i) {
    adj[cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
    adj[cd[el->edges[i].t] + d[el->edges[i].t]++] = el->edges[i].s;
  }

  BinaryHeap *heap = BH_MakeHeap(n, d);

  el->rank = (unsigned *)malloc(n * sizeof(unsigned));
  unsigned r = 0;
  for (unsigned i = 0; i < n; ++i) {
    KeyValuePair kv = BH_PopMin(heap);
    el->rank[kv.key] = n - (++r);
    for (unsigned j = cd[kv.key]; j < cd[kv.key + 1]; ++j) {
      BH_Update(heap, adj[j]);
    }
  }
  BH_FreeHeap(heap);
  free(d);
  free(cd);
  free(adj);
}

// Relabel each node to be its rank in degeneracy ordering.
// Each edge is oriented from the endpoint with lower (larger) rank to that with higher (smaller) rank.
void Relabel(EdgeList *el) {
  unsigned source, target, tmp;
  el->n = 0;
  for (unsigned i = 0; i < el->e; ++i) {
    source = el->rank[el->edges[i].s];
    target = el->rank[el->edges[i].t];
    if (source < target){
      tmp = source;
      source = target;
      target = tmp;
    }
    if (source + 1 > el->n){
      el->n = source + 1;
    }
    el->edges[i].s = source;
    el->edges[i].t = target;
  }
}

static int EdgeCmp(const void *a, const void *b) {
  if (((Edge *)a)->s != ((Edge *)b)->s) return ((Edge *)a)->s - ((Edge *)b)->s;
  return ((Edge *)a)->t - ((Edge *)b)->t;
}

// Building the directed graph
Graph *MakeGraph(EdgeList *el) {
  unsigned max;
  Graph *g = (Graph *)malloc(sizeof(Graph));
  unsigned *d = (unsigned *)calloc(el->n, sizeof(unsigned));

  qsort(el->edges, el->e, sizeof(Edge), EdgeCmp); // Sort the edges to achieve "reverse-order" enumeration. This sorting will ensure that all adjacency lists constructed for future subgraphs are sorted.

  for (unsigned i = 0; i < el->e; ++i) {
    ++d[el->edges[i].s];
  }

  g->cd = (unsigned *)malloc((el->n + 1) * sizeof(unsigned));
  g->cd[0] = 0;
  max = 0;
  for (unsigned i = 1; i < el->n + 1; ++i) {
    g->cd[i] = g->cd[i - 1] + d[i - 1];
    max = (max > d[i - 1]) ? max : d[i - 1];
    d[i - 1] = 0;
  }
  printf("Core value (max truncated degree) = %u\n", max);

  g->adj = (unsigned *)malloc(el->e * sizeof(unsigned));
  for (unsigned i = 0; i < el->e; ++i) {
    // printf("%u %u\n", el->edges[i].s, el->edges[i].t);
    g->adj[g->cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
  }

  free(d);
  g->core = max;
  g->n = el->n;
  free(el->rank);
  g->edges = el->edges;
  g->e = el->e;
  free(el);

  return g;
}

unsigned CountEdges(Graph *g, unsigned num_nodes, unsigned nodes[]) {
  unsigned m = 0;
  char *is_included = (char *)calloc(g->n, sizeof(char));
  for (unsigned i = 0; i < num_nodes; ++i)
    is_included[nodes[i]] = 1;
  for (unsigned i = 0; i < g->e; ++i)
    m += (is_included[g->edges[i].s] && is_included[g->edges[i].t]);
  return m;
}
