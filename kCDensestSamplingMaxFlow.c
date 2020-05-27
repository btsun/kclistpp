/*
Info:
  This program corresponds to the competitor "MaxFlow-Sampling" in the PVLDB
2020 paper.
  Feel free to use these lines as you wish.
  This program randomly selects a small fraction of k-cliques to store in main
memory, and uses a series of max-flow computations to find the densest subgraph
in the sampled hypergraph via a binary search.

To compile:
  "g++ kCDensestSamplingMaxFlow.c BinaryHeap.c Graph.c MaxFlowDouble.cpp -O3 -o kCDensestSamplingMaxFlow -lm -fopenmp"

To execute:
  "./kCDensestSamplingMaxFlow p k num_of_cliques_to_sample edgeListFileName".
  p is the number of threads.
  k is the size of a clique considered as in "k-clique".
  num_of_cliques_to_sample is the approximate number of k-cliques to sample. The
sampling probability will be set as this number divided by the total number of
k-cliques.
  edgeListFileName is the name of the file that contains the graph. Each line of
the file contains one edge represented by two integers separated by a space.

Output:
  The program will print relevant information throughout the computation,
including
  - the lower and upper bound obtained in each iteration of the binary
search;
  - the number of nodes of the densest subset in the sampled hypergraph;
  - the number of k-cliques and the k-clique density of that subset in the
original graph;
  - the number of max-flow computations;
  - the running time.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include "Graph.h"
#include "MaxFlowDouble.hpp"

static int UnsignedCmp(const void *a, const void *b) {
  return (long long)*(unsigned *)a - (long long)*(unsigned *)b;
}

inline int LargeRand() {
  if (RAND_MAX == 0x7fff)
    return (rand() << 15) | rand();
  return rand();
}

inline int GetRandMax() {
  if (RAND_MAX == 0x7fff)
    return 0x3fffffff;
  return RAND_MAX;
}

typedef enum {COUNT = 1, SAMPLING = 2, COUNT_IN_SUBGRAPH = 3} task_t;
static unsigned *original_graph_id_sg2g = NULL, *original_graph_id_g2sg = NULL; // to improve (???)
#pragma omp threadprivate(original_graph_id_g2sg, original_graph_id_sg2g)
unsigned *densest_subgraph_id_sg2g = NULL, *densest_subgraph_id_g2sg = NULL;

Subgraph *AllocSubgraph(Graph *g, unsigned char k) {
  Subgraph *sg = (Subgraph *)malloc(sizeof(Subgraph));
  sg->n = (unsigned *)calloc(k, sizeof(unsigned));
  sg->d = (unsigned **)malloc(k * sizeof(unsigned *));
  sg->adj = (unsigned *)malloc(g->core * g->core * sizeof(unsigned));
  sg->label = (unsigned char *)calloc(g->core, sizeof(unsigned char));
  sg->nodes = (unsigned **)malloc(k * sizeof(unsigned *));
  sg->core = g->core;
  for (unsigned i = 1; i < k; ++i){
    sg->d[i] = (unsigned *)malloc(g->core * sizeof(unsigned));
    sg->nodes[i] = (unsigned *)malloc(g->core * sizeof(unsigned));
  }
  return sg;
}

void MakeSubgraph(Graph *g, unsigned u, unsigned v, Subgraph *sg, unsigned char k, unsigned *id_sg2g, unsigned *id_g2sg, task_t task) {
  if (id_sg2g == NULL){
    id_g2sg = (unsigned *)malloc(g->n * sizeof(unsigned));
    id_sg2g = (unsigned *)malloc(g->core * sizeof(unsigned));
    for (unsigned i = 0; i < g->n; ++i) {
      id_g2sg[i] = UINT_MAX;
    }
  }

  for (unsigned i = 0; i < sg->n[k - 1]; ++i) {
    sg->label[i] = 0;
  }

  for (unsigned i = g->cd[v]; i < g->cd[v + 1]; ++i) { // For each out-neighbor of v
    id_g2sg[g->adj[i]] = UINT_MAX - 1;
  }

  unsigned j = 0;
  for (unsigned i = g->cd[u]; i < g->cd[u + 1]; ++i) { // For each out-neighbor of u
    unsigned x = g->adj[i];
    if (id_g2sg[x] == UINT_MAX - 1) {
      id_g2sg[x] = j;
      id_sg2g[j] = x;
      sg->label[j] = k - 2;
      sg->nodes[k - 2][j] = j;
      sg->d[k - 2][j] = 0; // New degrees
      ++j;
    }
  }
  sg->n[k - 2] = j;

  for (unsigned i = 0; i < sg->n[k - 2]; ++i) { // Reorder adjacency list and compute new degrees
    unsigned x = id_sg2g[i];
    for (unsigned l = g->cd[x]; l < g->cd[x + 1]; ++l) {
      unsigned y = g->adj[l];
      j = id_g2sg[y];
      if (j < UINT_MAX - 1) {
        sg->adj[sg->core * i + sg->d[k - 2][i]++] = j;
      }
    }
  }

  for (unsigned i = g->cd[v]; i < g->cd[v + 1]; ++i) {
    id_g2sg[g->adj[i]] = -1;
  }

  if (task == COUNT || task == SAMPLING) {
    original_graph_id_g2sg = id_g2sg;
    original_graph_id_sg2g = id_sg2g;
  } else {
    densest_subgraph_id_g2sg = id_g2sg;
    densest_subgraph_id_sg2g = id_sg2g;
  }
}

// ==========
// kCList: the clique-listing procedure
// ==========

unsigned num_of_cliques_to_sample;
unsigned sampled_cliques_reserved_size; // Maximum number of cliques for memory allocation; will increase if needed
unsigned *cknodes; // Nodes of a clique being formed
unsigned *ck; // List of all sampled cliques
unsigned *p_ckend; // Pointer to the end of ck[]
unsigned long long cnt_clique; // Number of cliques
unsigned long long cnt_sampled_clique; // Number of sampled cliques
double sampling_prob; // Sampling probability
unsigned *is_in_densest_subgraph; // Whether each node is in the densest subgraph
unsigned long long cnt_clique_in_densest_subgraph; // Number of cliques in the densest subgraph (without sampling)

bool clique_is_in_densest_subgraph(unsigned char clique_size) {
  for (unsigned i = 0; i < clique_size; ++i)
    if (!is_in_densest_subgraph[cknodes[i]])
      return false;
  return true;
}

void KCLIST_CliqueEnumThread(Subgraph *sg, unsigned char clique_size, unsigned char l, task_t task) {
  if (clique_size == 3) {
    for (unsigned i = 0; i < sg->n[1]; ++i) {
      unsigned u = sg->nodes[1][i];
      if (task == SAMPLING)
        cknodes[0] = original_graph_id_sg2g[u]; // When task == COUNT_IN_SUBGRAPH, cknodes is useless
      switch (task) {
        case COUNT: {
          ++cnt_clique;
          break;
        }
        case SAMPLING: {
          if (LargeRand() >= (GetRandMax() + 1LL) * sampling_prob) // Store this clique with probability sampling_prob
            break;
          #pragma omp critical
          {
            if (cnt_sampled_clique >= sampled_cliques_reserved_size) {
              sampled_cliques_reserved_size *= 2;
              ck = (unsigned *)realloc(ck, sampled_cliques_reserved_size * clique_size * sizeof(unsigned));
              p_ckend = ck + cnt_sampled_clique * clique_size;
            }
            for (unsigned j = 0; j < clique_size; ++j)
              *(p_ckend++) = cknodes[j];
            ++cnt_sampled_clique;
          }
          break;
        }
        case COUNT_IN_SUBGRAPH: {
          ++cnt_clique_in_densest_subgraph;
          break;
        }
      }
    }
    return;
  }
  if (l == 2) {
    for (unsigned i = 0; i < sg->n[2]; ++i) {
      unsigned u = sg->nodes[2][i];
      if (task == SAMPLING)
        cknodes[1] = original_graph_id_sg2g[u];
      for (unsigned j = u * sg->core, end = u * sg->core + sg->d[2][u]; j < end; ++j) {
        unsigned v = sg->adj[j];
        if (task == SAMPLING)
          cknodes[0] = original_graph_id_sg2g[v];
        switch (task) {
          case COUNT: {
            ++cnt_clique;
            break;
          }
          case SAMPLING: {
            if (LargeRand() > (GetRandMax() + 1LL) * sampling_prob) // Store this clique with probability sampling_prob
              break;
            #pragma omp critical
            {
              if (cnt_sampled_clique >= sampled_cliques_reserved_size) {
                sampled_cliques_reserved_size *= 2;
                ck = (unsigned *)realloc(ck, sampled_cliques_reserved_size * clique_size * sizeof(unsigned));
                p_ckend = ck + cnt_sampled_clique * clique_size;
              }
              for (unsigned k = 0; k < clique_size; ++k)
                *(p_ckend++) = cknodes[k];
              ++cnt_sampled_clique;
            }
            break;
          }
          case COUNT_IN_SUBGRAPH: {
            ++cnt_clique_in_densest_subgraph;
            break;
          }
        }
      }
    }
    return;
  }

  for (unsigned i = 0; i < sg->n[l]; ++i) { // Enumerate in reverse order. Very confusing! "++i" is actually the reverse order.
    unsigned u = sg->nodes[l][i];
    cknodes[l - 1] = original_graph_id_sg2g[u];

    sg->n[l - 1] = 0;
    unsigned end = u * sg->core + sg->d[l][u];
    for (unsigned j = u * sg->core; j < end; ++j) { // Relabel nodes and forming U'.
      unsigned v = sg->adj[j];
      if (sg->label[v] == l) {
        sg->label[v] = l - 1;
        sg->nodes[l - 1][sg->n[l - 1]++] = v;
        sg->d[l - 1][v] = 0; // New degrees
      }
    }
    for (unsigned j = 0; j < sg->n[l - 1]; ++j) { // Reorder adjacency list and compute new degrees
      unsigned v = sg->nodes[l - 1][j];
      for (unsigned k = sg->core * v, end = sg->core * v + sg->d[l][v]; k < end; ++k) {
        unsigned w = sg->adj[k];
        if (sg->label[w] == l - 1) {
          ++sg->d[l - 1][v];
        }
        else{
          sg->adj[k--] = sg->adj[--end];
          sg->adj[end] = w;
        }
      }
      qsort(sg->adj + sg->core * v, sg->d[l - 1][v], sizeof(unsigned), UnsignedCmp); // Sort the nodes in reverse order
    }

    KCLIST_CliqueEnumThread(sg, clique_size, l - 1, task);

    for (unsigned j = 0; j < sg->n[l - 1]; ++j) { // Restore labels
      unsigned v = sg->nodes[l - 1][j];
      sg->label[v] = l;
    }
  }
}

void KCLIST_CliqueEnum(Graph *g, unsigned char k, task_t task) {
  Subgraph *sg;
  switch (task) {
    case COUNT: {
      cnt_clique = 0;
      break;
    }
    case SAMPLING: {
      cnt_sampled_clique = 0;
      sampled_cliques_reserved_size = 1.1 * num_of_cliques_to_sample;
      sampling_prob = (num_of_cliques_to_sample < cnt_clique) ? (double)num_of_cliques_to_sample / cnt_clique : 1;
      p_ckend = ck = (unsigned *)malloc(sampled_cliques_reserved_size * k * sizeof(unsigned));
      break;
    }
    case COUNT_IN_SUBGRAPH: {
      cnt_clique_in_densest_subgraph = 0;
      break;
    }
  }
  if (task == COUNT || task == SAMPLING) {
    #pragma omp parallel private(sg) reduction(+: cnt_sampled_clique)
    {
      cknodes = (unsigned *)malloc(k * sizeof(unsigned));
      sg = AllocSubgraph(g, k);
      #pragma omp for schedule(dynamic, 1) nowait
      for(unsigned i = 0; i < g->e; ++i) {
        cknodes[k - 1] = g->edges[i].s;
        cknodes[k - 2] = g->edges[i].t;
        MakeSubgraph(g, g->edges[i].s, g->edges[i].t, sg, k, original_graph_id_sg2g, original_graph_id_g2sg, task);
        KCLIST_CliqueEnumThread(sg, k, k - 2, task);
      }
      FreeSubgraph(sg, k);
    }
  } else {
    cknodes = (unsigned *)malloc(k * sizeof(unsigned));
    sg = AllocSubgraph(g, k);
    densest_subgraph_id_g2sg = densest_subgraph_id_sg2g = NULL;
    for (unsigned i = 0; i < g->e; ++i) {
      MakeSubgraph(g, g->edges[i].s, g->edges[i].t, sg, k, densest_subgraph_id_sg2g, densest_subgraph_id_g2sg, task);
      KCLIST_CliqueEnumThread(sg, k, k - 2, task);
    }
    free(densest_subgraph_id_g2sg);
    free(densest_subgraph_id_sg2g);
    free(cknodes);
    FreeSubgraph(sg, k);
  }
  switch (task) {
    case COUNT: {
      printf("Number of %u-cliques: %llu\n", k, cnt_clique);
      break;
    }
    case SAMPLING: {
      ck = (unsigned *)realloc(ck, cnt_sampled_clique * k * sizeof(unsigned));
      printf("Number of sampled %u-cliques: %llu\n", k, cnt_sampled_clique);
      break;
    }
    case COUNT_IN_SUBGRAPH: {
      break;
    }
  }
}

EdgeList *MakeDensestSubgraphEdgeList(Graph *g, const unsigned char k, const unsigned densest_subset_size) {
  EdgeList *el = (EdgeList *)malloc(sizeof(EdgeList));
  unsigned *new_id = (unsigned *)malloc(g->n * sizeof(unsigned));
  for (unsigned i = 0, j = 0; i < g->n; ++i) {
    if (is_in_densest_subgraph[i]) {
      new_id[i] = j++;
    }
  }
  el->n = densest_subset_size;
  el->e = 0;
  for (unsigned i = 0; i < g->e; ++i)
    el->e += (is_in_densest_subgraph[g->edges[i].s] && is_in_densest_subgraph[g->edges[i].t]);
  el->edges = (Edge *)malloc(el->e * sizeof(Edge));
  for (unsigned i = 0, j = 0; i < g->e; ++i) {
    if (is_in_densest_subgraph[g->edges[i].s] && is_in_densest_subgraph[g->edges[i].t]) {
      el->edges[j].s = new_id[g->edges[i].s];
      el->edges[j].t = new_id[g->edges[i].t];
      ++j;
    }
  }
  free(new_id);
  return el;
}

void Solve(const unsigned char k, Graph *g, time_t *p_t1) {
  const double EPS = 1e-9;

  // Count the number of clqiues
  KCLIST_CliqueEnum(g, k, COUNT);
  time_t t1 = *p_t1;
  time_t t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  *p_t1 = t2;

  // Sample k-cliques
  KCLIST_CliqueEnum(g, k, SAMPLING);
  t1 = *p_t1;
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  *p_t1 = t2;

  // Compute k-clique degree in the sampled hypergraph
  unsigned *ckdeg = (unsigned *)calloc(g->n, sizeof(unsigned));
  is_in_densest_subgraph = (unsigned *)calloc(g->n, sizeof(unsigned));
  for (unsigned long long i = 0; i < cnt_sampled_clique * k; ++i)
    ++ckdeg[ck[i]];

  // Construct the max-flow network
  double l = (double)cnt_sampled_clique / g->n;
  double u = 1;
  for (unsigned i = 1; i < k; ++i)
    u *= (double)(g->n - i) / (i + 1);
  Network network;
  Network::Vertex s = network.AddVertex();
  Network::Vertex t = network.AddVertex();
  vector<Network::Vertex> L;
  vector<Network::Edge> edgesToCheck; // Examine these edges later to determine the minimum cut

  // Construct the network as described in Algorithm 6 of the WWW 2015 paper by Trourakakis
  for (unsigned i = 0; i < g->n; ++i) {
    L.push_back(network.AddVertex());
    network.AddEdge(s, L[i], ckdeg[i]);
    // The capacity of the following edge will be modified for each binary search step, and the residual capacity will be checked after max-flow computation
    edgesToCheck.push_back(network.AddEdge(L[i], t, 0));
  }
  for (unsigned long long i = 0; i < cnt_sampled_clique; ++i) {
    Network::Vertex v = network.AddVertex();
    for (unsigned j = 0; j < k; ++j) {
      network.AddEdge(L[ck[i * k + j]], v, 1);
      network.AddEdge(v, L[ck[i * k + j]], k - 1);
    }
  }
  /*
  // Construct the network as described in Construction B of the KDD 2015 paper by Mitzenmacher et al.
  for (unsigned i = 0; i < g->n; ++i) {
    L.push_back(network.AddVertex());
    // The capacity of the following edge will be modified for each binary search step, and the residual capacity will be checked after max-flow computation
    edgesToCheck.push_back(network.AddEdge(s, L[i], 0));
  }
  for (unsigned long long i = 0; i < cnt_sampled_clique; ++i) {
    Network::Vertex v = network.AddVertex();
    for (unsigned j = 0; j < k; ++j) {
      network.AddEdge(L[ck[i * k + j]], v, u); // Infinite capacity
    }
    network.AddEdge(v, t, 1);
  }
  */
  // Binary search
  unsigned cnt_max_flow = 0;
  while (u - l >= 1.0 / g->n / (g->n - 1) && cnt_max_flow < 200) {
    printf("[Lower Bound, Upper Bound] = [%.12f, %.12f]\n", l, u);
    // fprintf(ofp, "%.12f\t%.12f\t%u\t%ld\n", l, u, cnt_max_flow, time(NULL) - t0);
    // fflush(ofp);
    double guess = (l + u) / 2;
    for (unsigned i = 0; i < g->n; ++i) {
      network.capacity[edgesToCheck[i]] = k * guess; // This line is for Algorithm 6 of the WWW 2015 paper by Trourakakis
      // network.capacity[edgesToCheck[i]] = guess; // This line is for Construction B of the KDD 2015 paper by Mitzenmacher et al.
    }

    double max_flow = network.MaxFlow(s, t);
    if (max_flow < cnt_sampled_clique * k * (1 - EPS)) { // This line is for Algorithm 6 of the WWW 2015 paper by Trourakakis
    // if (max_flow < cnt_sampled_clique * (1 - EPS)) { // This line is for Construction B of the KDD 2015 paper by Mitzenmacher et al.
      l = guess;
    }
    else
      u = guess;
    ++cnt_max_flow;
  }
  printf("Binary search completed.\n");
  t1 = *p_t1;
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  *p_t1 = t2;

  // Identify the densest subgraph
  unsigned cnt_node_in_densest_subgraph = 0;
  for (unsigned i = 0; i < g->n; ++i) {
    network.capacity[edgesToCheck[i]] = k * l; // This line is for Algorithm 6 of the WWW 2015 paper by Trourakakis
    // network.capacity[edgesToCheck[i]] = l; // This line is for Construction B of the KDD 2015 paper by Mitzenmacher et al.
  }
  network.MaxFlow(s, t);
  network.BfsOnResidualGraph(s);
  for (unsigned i = 0; i < g->n; ++i) {
    is_in_densest_subgraph[i] = network.color[L[i]] == boost::black_color; // This line is for Algorithm 6 of the WWW 2015 paper by Trourakakis
    // is_in_densest_subgraph[i] = network.color[L[i]] == boost::white_color; // This line is for Construction B of the KDD 2015 paper by Mitzenmacher et al.
    cnt_node_in_densest_subgraph += is_in_densest_subgraph[i];
  }

  // Count the number of k-cliques in the subgraph of the original graph induced by the densest subset
  EdgeList *el = MakeDensestSubgraphEdgeList(g, k, cnt_node_in_densest_subgraph);
  SortByCore(el);
  Relabel(el);
  Graph *p_densest_subgraph = MakeGraph(el);
  KCLIST_CliqueEnum(p_densest_subgraph, k, COUNT_IN_SUBGRAPH);

  printf("The densest subgraph has %d nodes and %lld %d-cliques.\n", cnt_node_in_densest_subgraph, cnt_clique_in_densest_subgraph, k);
  printf("Clique density is %.9f\n", (double)cnt_clique_in_densest_subgraph / cnt_node_in_densest_subgraph);
  printf("Number of max-flows = %d.\n", cnt_max_flow);
}

int main(int argc, char **argv) {
  EdgeList *el;
  Graph *g;
  unsigned num_threads = atoi(argv[1]);
  unsigned char k = atoi(argv[2]);
  num_of_cliques_to_sample = atoi(argv[3]);
  char *file_name = argv[4];
  omp_set_num_threads(num_threads);

  time_t t0, t1, t2;
  t0 = t1 = time(NULL);

  printf("Reading edgelist from file %s\n", file_name);
  el = ReadEdgeList(file_name);
  printf("Number of nodes = %u\n", el->n);
  printf("Number of edges = %u\n", el->e);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n",(t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  printf("Building the graph structure\n");
  SortByCore(el); // Do core decomposition and render degeneracy ordering to the nodes
  Relabel(el);
  g = MakeGraph(el);
  printf("Number of nodes (degree > 0) = %u\n", g->n);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  Solve(k, g, &t1);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  FreeGraph(g);
  printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));
  //fprintf(ofp, "%ld\n", t2 - t0);
  //fclose(ofp);
  return 0;
}
