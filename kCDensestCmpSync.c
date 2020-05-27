/*
Info:
  This program corresponds to "Sim-kClist++" in the PVLDB 2020 paper.
  Feel free to use these lines as you wish.
  This program iterates over all k-cliques for many rounds to run the Frank-
Wolfe based algortihm using linear memory and report approximate maximum
k-clique density.
  This program can handle both the case k = 2 and the case k >= 3.

To compile:
  "gcc kCDensestCmpSync.c BinaryHeap.c Graph.c -O3 -o kCDensestCmpSync -lm -fopenmp"

To execute:
  "./kCDensest p T k edgeListFileName tag"
  p is the number of threads.
  T is the number of iterations of the Frank-Wolfe algorithm (will be rounded
down to the nearest power of 2).
  k is the size of a clique considered as in "k-clique".
  edgeListFileName is the name of the file that contains the graph. Each line of
the file contains one edge represented by two integers separated by a space.
  tag is a string specifying the dataset (e.g., "dblp"), which is used to
generate the output file name.

Output:
  Evolution of the approximate k-clique densest subgraph. One record per line,
containing
  - the number of iterations completed
  - the number of nodes in the approximate k-clique densest subgraph;
  - the number of edges in the approximate k-clique densest subgraph;
  - the k-clique density of the approximate k-clique densest subgraph;
  - the computed upper bound on the maximum k-clique density;
  - the time elapsed since the beginning of the execution.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include "Graph.h"

static int UnsignedCmp(const void *a, const void *b) {
  return (long long)*(unsigned *)a - (long long)*(unsigned *)b;
}

static int EdgeCmp(const void *a, const void *b) {
  if (((Edge *)a)->s != ((Edge *)b)->s) return ((Edge *)a)->s - ((Edge *)b)->s;
  return ((Edge *)a)->t - ((Edge *)b)->t;
}

Subgraph *AllocSubgraph(Graph *g, unsigned char k) {
  Subgraph *sg = malloc(sizeof(Subgraph));
  sg->n = calloc(k, sizeof(unsigned));
  sg->d = malloc(k * sizeof(unsigned *));
  sg->adj = malloc(g->core * g->core * sizeof(unsigned));
  sg->label = calloc(g->core, sizeof(unsigned char));
  sg->nodes = malloc(k * sizeof(unsigned *));
  sg->core = g->core;
  for (unsigned i = 0; i < k; ++i){
    sg->d[i] = malloc(g->core * sizeof(unsigned));
    sg->nodes[i] = malloc(g->core * sizeof(unsigned));
  }
  return sg;
}

static unsigned *id_sg2g = NULL, *id_g2sg = NULL; // to improve (???)
#pragma omp threadprivate(id_g2sg, id_sg2g)

void MakeSubgraph(Graph *g, Edge edge, Subgraph *sg, unsigned char k) {
  unsigned u = edge.s, v = edge.t;

  if (id_sg2g == NULL){
    id_g2sg = malloc(g->n * sizeof(unsigned));
    id_sg2g = malloc(g->core * sizeof(unsigned));
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
double *rho;
double *rho_add;
double *rho_tentative;
unsigned *level;
unsigned *reordered;
typedef enum {FRANK_WOLFE_INIT = 1, FRANK_WOLFE = 2, PAVA_PREPROCESS = 3} task_t;

void AllocCdf(Graph *g, unsigned k) {
  rho = malloc(g->n * sizeof(double));
  rho_add = malloc(g->n * sizeof(double));
  rho_tentative = malloc(g->n * sizeof(double));
  level = malloc(g->n * sizeof(unsigned));
  reordered = malloc(g->n * sizeof(unsigned));
}

void CDF_FrankWolfeUpdateRates(int clique_size) {
  unsigned node_getting_weight = cknodes[0];
  for (unsigned i = 1; i < clique_size; ++i) {
    if (rho[node_getting_weight] > rho[cknodes[i]])
      node_getting_weight = cknodes[i];
  }
  #pragma omp atomic
  rho_add[node_getting_weight] += 1.0;
}

void CDF_CliqueEnumThread(Subgraph *sg,
                          unsigned char clique_size,
                          unsigned char l,
                          task_t task,
                          unsigned node_getting_weight) {
  switch (task) {
    case FRANK_WOLFE_INIT: {
      if (clique_size == 2) {
        #pragma omp atomic
        rho[cknodes[0]] += 0.5;
        #pragma omp atomic
        rho[cknodes[1]] += 0.5;
        return;
      }
      if (clique_size == 3) {
        // Set initial rho values: each member of each clique gets 1 / k
        rho[cknodes[2]] += (double)sg->n[1] / (double)clique_size;
        rho[cknodes[1]] += (double)sg->n[1] / (double)clique_size;
        for (unsigned i = 0; i < sg->n[1]; ++i) {
          unsigned u = sg->nodes[1][i];
          rho[id_sg2g[u]] += 1.0 / 3.0;
        }
        return;
      }
      if (l == 2) {
        // Set initial rho values: each member of each clique gets 1 / k
        for (unsigned i = 0; i < sg->n[2]; ++i) {
          unsigned u = sg->nodes[2][i];
          for (unsigned j = u * sg->core, end = u * sg->core + sg->d[2][u]; j < end; ++j) {
            unsigned v = sg->adj[j];
            rho[id_sg2g[v]] += 1.0 / clique_size;
          }
          rho[id_sg2g[u]] += (double)sg->d[2][u] / (double)clique_size;
          for (unsigned j = 2; j < clique_size; ++j)
            rho[cknodes[j]] += (double)sg->d[2][u] / (double)clique_size;
        }
        return;
      }
      break;
    }
    case FRANK_WOLFE: {
      if (clique_size == 2) {
        CDF_FrankWolfeUpdateRates(clique_size);
        return;
      }
      if (clique_size == 3) {
        // Adjust rho values: pick the node with least rho value from each clique
        for (unsigned i = 0; i < sg->n[1]; ++i) {
          unsigned u = sg->nodes[1][i];
          cknodes[0] = id_sg2g[u];
          CDF_FrankWolfeUpdateRates(clique_size);
        }
        return;
      }
      if (l == 2) {
        // Adjust rho values: pick the node with least rho value from each clique
        for(unsigned i = 0; i < sg->n[2]; ++i) { // List all edges
          unsigned u = sg->nodes[2][i];
          cknodes[1] = id_sg2g[u];
          for (unsigned j = u * sg->core, end = u * sg->core + sg->d[2][u]; j < end; ++j) {
            unsigned v = sg->adj[j];
            cknodes[0] = id_sg2g[v];
            CDF_FrankWolfeUpdateRates(clique_size);
          }
        }
        return;
      }
      break;
    }
    case PAVA_PREPROCESS: {
      if (clique_size == 2) {
        #pragma omp atomic
        rho_tentative[level[node_getting_weight]] += 1.0;
        return;
      }
      if (clique_size == 3) {
        for (unsigned i = 0; i < sg->n[1]; ++i) {
          unsigned u = sg->nodes[1][i];
          unsigned node_getting_weight_2 = level[id_sg2g[u]] > level[node_getting_weight] ? id_sg2g[u] : node_getting_weight;
          #pragma omp atomic
          rho_tentative[level[node_getting_weight_2]] += 1.0;
        }
        return;
      }
      if (l == 2) {
        for (unsigned i = 0; i < sg->n[2]; ++i) { // List all edges
          unsigned u = sg->nodes[2][i];
          unsigned node_getting_weight_2 = level[id_sg2g[u]] > level[node_getting_weight] ? id_sg2g[u] : node_getting_weight;
          for (unsigned j = u * sg->core, end = u * sg->core + sg->d[2][u]; j < end; ++j) {
            unsigned v = sg->adj[j];
            unsigned node_getting_weight_3 = level[id_sg2g[v]] > level[node_getting_weight_2] ? id_sg2g[v] : node_getting_weight_2;
            #pragma omp atomic
            rho_tentative[level[node_getting_weight_3]] += 1.0;
          }
        }
        return;
      }
      break;
    }
  }

  for(unsigned i = 0; i < sg->n[l]; ++i) { // Enumerate in reverse order. Very confusing! "++i" is actually the reverse order.
    unsigned u = sg->nodes[l][i];
    cknodes[l - 1] = id_sg2g[u];
    unsigned node_getting_weight_2;
    switch (task) {
      case FRANK_WOLFE: {
        node_getting_weight_2 = rho[id_sg2g[u]] < rho[node_getting_weight] ? id_sg2g[u] : node_getting_weight;
        break;
      }
      case PAVA_PREPROCESS: {
        node_getting_weight_2 = level[id_sg2g[u]] > level[node_getting_weight] ? id_sg2g[u] : node_getting_weight;
        break;
      }
    }

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

    CDF_CliqueEnumThread(sg, clique_size, l - 1, task, node_getting_weight_2);

    for (unsigned j = 0; j < sg->n[l - 1]; ++j) { // Restore labels
      unsigned v = sg->nodes[l - 1][j];
      sg->label[v] = l;
    }
  }
}

void CDF_CliqueEnum(Graph *g, unsigned char k, task_t task) {
  Subgraph *sg;
  #pragma omp parallel private(sg)
  {
    cknodes = malloc(k * sizeof(unsigned));
    sg = AllocSubgraph(g, k);
    #pragma omp for schedule(dynamic, 1) nowait
    for(unsigned i = 0; i < g->e; ++i) {
      cknodes[k - 1] = g->edges[i].s;
      cknodes[k - 2] = g->edges[i].t;
      switch (task) {
        case FRANK_WOLFE: {
          unsigned node_getting_weight = rho[cknodes[k - 2]] < rho[cknodes[k - 1]] ? cknodes[k - 2] : cknodes[k - 1];
          MakeSubgraph(g, g->edges[i], sg, k);
          CDF_CliqueEnumThread(sg, k, k - 2, FRANK_WOLFE, node_getting_weight);
          break;
        }
        case PAVA_PREPROCESS: {
          unsigned node_getting_weight = level[cknodes[k - 2]] > level[cknodes[k - 1]] ? cknodes[k - 2] : cknodes[k - 1];
          MakeSubgraph(g, g->edges[i], sg, k);
          CDF_CliqueEnumThread(sg, k, k - 2, PAVA_PREPROCESS, node_getting_weight);
          break;
        }
      }
    }
    FreeSubgraph(sg, k);
  }
}

typedef struct {
  unsigned n; // Number of nodes
  unsigned m; // Number of edges
  double density;
  double ub; // An upper bound of maximum density
} DensestSubsetInfo;

static int CDF_NodeCmp(const void *a, const void *b) {
  unsigned long long x = rho[*(const unsigned *)a];
  unsigned long long y = rho[*(const unsigned *)b];
  if (x > y) return -1;
  if (x < y) return 1;
  return 0;
}

DensestSubsetInfo CDF_FindDensestSubset(Graph *g, unsigned char k) {
  DensestSubsetInfo info;
  info.density = -1;
  for (unsigned i = 0; i < g->n; ++i)
    reordered[i] = i;
  qsort(reordered, g->n, sizeof(unsigned), CDF_NodeCmp); // Reorder the nodes by decreasing rho values
  for (unsigned i = 0; i < g->n; ++i)
    level[reordered[i]] = i;
//  for (int i = 0; i < g->n; ++i)
//    printf("level[%u] = %u\n", i, level[i]);
  CDF_CliqueEnum(g, k, PAVA_PREPROCESS);
  // Find the approximate maximum density
  double sum = 0;
  for (unsigned i = 0; i < g->n; ++i) {
    sum += rho_tentative[i];
    if (sum / (i + 1) > info.density) {
      info.density = (double)sum / (i + 1);
      info.n = i + 1;
    }
  }
  // Count the number of edges
  info.m = CountEdges(g, info.n, reordered);
  // Compute an upper bound of maximum density
  sum = 0;
  info.ub = 0;
  double ip1ck = 1; // (i + 1) choose k
  for (unsigned i = 0; i < g->n; ++i) {
    sum += rho[reordered[i]];
    if (i + 1 == k)
      ip1ck = 1;
    else if (i + 1 > k)
      ip1ck = (ip1ck * (i + 1)) / (i + 1 - k);
    if (ip1ck < sum)
      info.ub = ip1ck / (i + 1);
    else {
      if (info.ub < sum / (i + 1))
        info.ub = sum / (i + 1);
      break;
    }
  }
  return info;
}

void CDF_Main(unsigned char k, Graph *g, unsigned num_iter, char *output_file_name, time_t t0) {
  FILE *ofp = fopen(output_file_name, "w");
  fprintf(ofp, "[Number of Iterations]\t[Number of Nodes]\t[Number of Edges]\t[k-Clique Density]\t[Upper Bound of Maximum Density]\t[Time (seconds)]\n");
  AllocCdf(g, k);
  qsort(g->edges, g->e, sizeof(Edge), EdgeCmp); // Sort the edges to achieve "reverse-order" enumeration (edge parallel)
  // for (unsigned u = 0; u < el->n; ++u) // Sort the edges to achieve "reverse-order" enumeration (node parallel)
  //   qsort(g->adj + g->cd[u], d[u], sizeof(unsigned), UnsignedCmp);
  for (unsigned i = 0; i < g->n; ++i)
    rho[i] = 0;
  for (unsigned T = 1, t = 1; T <= num_iter; T <<= 1) {
    // Step 1: run the Frank-Wolfe based algorithm for num_iter rounds
    for (; t <= T; ++t) {
      if (t % 10 == 0)
        fprintf(stderr, "Run round %u...\n", t);
      for (unsigned i = 0; i < g->n; ++i)
        rho_add[i] = 0;
      CDF_CliqueEnum(g, k, FRANK_WOLFE);
      for (unsigned i = 0; i < g->n; ++i)
        rho[i] = (t * rho[i] + 2.0 * rho_add[i]) / (t + 2.0);
    }
    // Step 2: give a tentative decomposition
    for (unsigned i = 0; i < g->n; ++i)
      rho_tentative[i] = 0;
    DensestSubsetInfo info = CDF_FindDensestSubset(g, k);
    time_t t1 = time(NULL);
    fprintf(stderr, "Approximate densest subgraph: %u nodes, density = %f, upper bound = %f. %ld seconds.\n", info.n, info.density, info.ub, t1 - t0);
    fprintf(ofp, "%u\t%u\t%u\t%.12f\t%.12f\t%ld\n", T, info.n, info.m, info.density, info.ub, t1 - t0);
    fflush(ofp);
  }
  fclose(ofp);
  /*FILE *ofp = fopen("rates.txt", "w");
  for (unsigned i = 0; i < g->n; ++i)
    fprintf(ofp, "r[%u] = %.12f\n", reordered[i], rho[reordered[i]]);
  fclose(ofp);*/
}

int main(int argc, char **argv) {
  EdgeList *el;
  Graph *g;
  unsigned num_iter = atoi(argv[2]);
  unsigned char k = atoi(argv[3]);
  char *file_name = argv[4];
  omp_set_num_threads(atoi(argv[1]));

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
  char output_file_name[100] = "stat_cmpsync_";
  strcat(output_file_name, argv[5]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[1]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[2]);
  strcat(output_file_name, "_");
  strcat(output_file_name, argv[3]);
  strcat(output_file_name, ".txt");
  CDF_Main(k, g, num_iter, output_file_name, t0);
  t2 = time(NULL);
  printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
  t1 = t2;

  FreeGraph(g);
  printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

  return 0;
}
