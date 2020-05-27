/*
Info:
  This program corresponds to the competitor "MaxFlow" in the PVLDB 2020 paper.
  Feel free to use these lines as you wish.
  This program implements the exact algorithm for the densest subgraph problem
and the k-clique densest subgraph problem. The algorithm for the former is
proposed by Goldberg in 1984, while the algorithm for the latter is proposed by
Trourakakis (WWW 2015). Actually the latter also applies to the densest subgraph
problem, i.e., the case k = 2. Both are based on binary search and max flow.

To compile:
  "g++ kCDensestG.cpp BinaryHeap.c Graph.c MaxFlowDouble.cpp -O3 -o kCDensestG -lm -fopenmp"

To execute:
  "./kCDensestG p k edgeListFileName tag".
  p is the number of threads.
  k is the size of a clique considered as in "k-clique".
  edgeListFileName is the name of the file that contains the graph. Each line of
the file contains one edge represented by two integers separated by a space.
  tag is a string specifying the dataset (e.g., "dblp"), which is used to
generate the output file name.

Output:
  The program will print the following for each iteration of the binary search:
  - the lower and upper bound on the maximum k-clique density;
  - the time elapsed since the beginning of the execution.
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

unsigned MAX_CLIQUES = 100000000; // Maximum number of cliques for memory allocation; will increase if needed

static int UnsignedCmp(const void *a, const void *b) {
  return (long long)*(unsigned *)a - (long long)*(unsigned *)b;
}

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

static unsigned *id_sg2g = NULL, *id_g2sg = NULL; // to improve (???)
#pragma omp threadprivate(id_g2sg, id_sg2g)

void MakeSubgraph(Graph *g, unsigned u, unsigned v, Subgraph *sg, unsigned char k) {
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
}

// Clique-density-friendly decomposition

unsigned *cknodes; // Nodes of a clique
#pragma omp threadprivate(cknodes)
unsigned *ck; // List of all cliques
unsigned *p_ckend; // Pointer to the end of ck[]
unsigned long long cnt_clique;

void CDF_CliqueEnumThread(Subgraph *sg,
                          unsigned char clique_size,
                          unsigned char l) {
  if (clique_size == 3) {
    for (unsigned i = 0; i < sg->n[1]; ++i) {
      unsigned u = sg->nodes[1][i];
      cknodes[0] = id_sg2g[u];
      #pragma omp critical
      {
        if (cnt_clique >= MAX_CLIQUES) {
          MAX_CLIQUES *= 2;
          ck = (unsigned *)realloc(ck, MAX_CLIQUES * clique_size * sizeof(unsigned));
          p_ckend = ck + cnt_clique * clique_size;
        }
        for (unsigned j = 0; j < clique_size; ++j)
          *(p_ckend++) = cknodes[j];
        ++cnt_clique;
      }
    }
    return;
  }
  if (l == 2) {
    for (unsigned i = 0; i < sg->n[2]; ++i) {
      unsigned u = sg->nodes[2][i];
      cknodes[1] = id_sg2g[u];
      for (unsigned j = u * sg->core, end = u * sg->core + sg->d[2][u]; j < end; ++j) {
        unsigned v = sg->adj[j];
        cknodes[0] = id_sg2g[v];
        #pragma omp critical
        {
          if (cnt_clique >= MAX_CLIQUES) {
            MAX_CLIQUES *= 2;
            ck = (unsigned *)realloc(ck, MAX_CLIQUES * clique_size * sizeof(unsigned));
            p_ckend = ck + cnt_clique * clique_size;
          }
          for (unsigned k = 0; k < clique_size; ++k)
            *(p_ckend++) = cknodes[k];
          ++cnt_clique;
        }
      }
    }
    return;
  }

  for (unsigned i = 0; i < sg->n[l]; ++i) { // Enumerate in reverse order. Very confusing! "++i" is actually the reverse order.
    unsigned u = sg->nodes[l][i];
    cknodes[l - 1] = id_sg2g[u];

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

    CDF_CliqueEnumThread(sg, clique_size, l - 1);

    for (unsigned j = 0; j < sg->n[l - 1]; ++j) { // Restore labels
      unsigned v = sg->nodes[l - 1][j];
      sg->label[v] = l;
    }
  }
}

void CDF_CliqueEnum(Graph *g, unsigned char k) {
  Subgraph *sg;
  cnt_clique = 0;
  p_ckend = ck = (unsigned *)malloc(MAX_CLIQUES * k * sizeof(unsigned));
  #pragma omp parallel private(sg) reduction(+: cnt_clique)
  {
    cknodes = (unsigned *)malloc(k * sizeof(unsigned));
    sg = AllocSubgraph(g, k);
    #pragma omp for schedule(dynamic, 1) nowait
    for(unsigned i = 0; i < g->e; ++i) {
      cknodes[k - 1] = g->edges[i].s;
      cknodes[k - 2] = g->edges[i].t;
      MakeSubgraph(g, g->edges[i].s, g->edges[i].t, sg, k);
      CDF_CliqueEnumThread(sg, k, k - 2);
    }
    FreeSubgraph(sg, k);
  }
  ck = (unsigned *)realloc(ck, cnt_clique * k * sizeof(unsigned));
}

void Solve_Main(const unsigned char k, Graph *g, FILE *ofp, time_t t0) {
  fprintf(ofp, "[Lower Bound]\t[Upper Bound]\t[Number of Max-Flows]\t[Time (seconds)]\n");
  if (k == 2) {
    ck = (unsigned *)malloc((unsigned long long)(g->e) * k * sizeof(unsigned));
    cnt_clique = g->e;
    for (unsigned long long i = 0; i < (unsigned long long)(g->e); ++i) {
      ck[i << 1] = g->edges[i].s;
      ck[(i << 1) + 1] = g->edges[i].t;
    }
  } else {
    CDF_CliqueEnum(g, k); // Collect all k-cliques
  }
  // Compute k-clique degree
  unsigned *ckdeg = (unsigned *)calloc(g->n, sizeof(unsigned));
  for (unsigned long long i = 0; i < cnt_clique * k; ++i)
    ++ckdeg[ck[i]];
  // Binary search
  double l = (double)cnt_clique / g->n;
  double u = 1;
  for (unsigned i = 1; i < k; ++i)
    u *= (double)(g->n - i) / (i + 1);
  if (u > cnt_clique)
    u = cnt_clique;
  unsigned cnt_max_flow = 0;
  while (u - l >= 1.0 / g->n / (g->n - 1)) {
    fprintf(stderr, "[Lower Bound, Upper Bound] = [%.12f, %.12f]\n", l, u);
    fprintf(ofp, "%.12f\t%.12f\t%u\t%ld\n", l, u, cnt_max_flow, time(NULL) - t0);
    fflush(ofp);
    double guess = (l + u) / 2;
    if (guess == l || guess == u) {
      fprintf(ofp, "Cannot achieve desired precision!!!\n");
      fprintf(stderr, "Cannot achieve desired precision!!!\n");
      break;
    }
    Network network;
    Network::Vertex s = network.AddVertex();
    Network::Vertex t = network.AddVertex();
    vector<Network::Vertex> L;
    if (k == 2) { // Construct the network as described in the paper "Finding a maximum density subgraph" by Goldberg in 1984
      for (unsigned i = 0; i < g->n; ++i) {
        L.push_back(network.AddVertex());
        network.AddEdge(s, L[i], g->e);
        network.AddEdge(L[i], t, g->e + 2 * guess - ckdeg[i]);
      }
      for (unsigned i = 0; i < g->e; ++i) {
        network.AddEdge(L[g->edges[i].s], L[g->edges[i].t], 1);
        network.AddEdge(L[g->edges[i].t], L[g->edges[i].s], 1);
      }
      double max_flow = network.MaxFlow(s, t);
      if (max_flow < (double)g->n * g->e * 0.999999)
        l = guess;
      else
        u = guess;
    } else { // Construct the network as described in Algorithm 6 of the WWW 2015 paper by Trourakakis
      for (unsigned i = 0; i < g->n; ++i) {
        L.push_back(network.AddVertex());
        network.AddEdge(s, L[i], ckdeg[i]);
        network.AddEdge(L[i], t, k * guess);
      }
      for (unsigned long long i = 0; i < cnt_clique; ++i) {
        Network::Vertex v = network.AddVertex();
        for (unsigned j = 0; j < k; ++j) {
          network.AddEdge(L[ck[i * k + j]], v, 1);
          network.AddEdge(v, L[ck[i * k + j]], k - 1);
        }
      }
      double max_flow = network.MaxFlow(s, t);
      if (max_flow < cnt_clique * k * 0.999999)
        l = guess;
      else
        u = guess;
    }
    ++cnt_max_flow;
  }
  fprintf(stderr, "Maximum density = %.12f. Number of max-flows = %u.\n", l, cnt_max_flow);
  fprintf(ofp, "%.12f\t%.12f\t%u\t%ld\n", l, u, cnt_max_flow, time(NULL) - t0);
}

int main(int argc, char **argv) {
  EdgeList *el;
  Graph *g;
  unsigned char k = atoi(argv[2]);
  char *file_name = argv[3];
  unsigned num_threads = atoi(argv[1]);
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

  printf("Iterate over all cliques\n");
  char output_file_name[100] = "stat_exacttheirs_";
  strcat(output_file_name, argv[4]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[1]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[2]);
  strcat(output_file_name, ".txt");
  FILE *ofp = fopen(output_file_name, "w");
  Solve_Main(k, g, ofp, t0);

  printf("Number of %u-cliques: %llu\n", k, cnt_clique);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  FreeGraph(g);
  printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));
  fprintf(ofp, "Total time (seconds) = %ld\n", t2 - t0);
  fclose(ofp);
  return 0;
}
