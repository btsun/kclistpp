/*
Info:
  This program corresponds to "Seq-Sampling++" in the PVLDB 2020 paper.
  Feel free to use these lines as you wish.
  This program iterates over all k-cliques, randomly saves a small part of them
in the main memory, iterates over these sampled k-cliques for many rounds and
report the approximate maximum k-clique density.
  Note that this program can only handle k >= 3, i.e., k = 2 is not supported.

To compile:
  "gcc kCDensestSamplingKclistpp.c BinaryHeap.c Graph.c -O3 -o kCDensestSamplingKclistpp -lm -fopenmp"

To execute:
  "./kCDensestSamplingKclistpp p T k edgeListFileName"
  p is the number of threads.
  T is the number of iterations of the "++" operation (will be rounded down to
the nearest power of 2).
  k is the size of a clique considered as in "k-clique". It must be at least 3.
  edgeListFileName is the name of the file that contains the graph. Each line of
the file contains one edge represented by two integers separated by a space.

Output:
  Evolution of the approximate k-clique densest subgraph. One record per line,
containing
  - the number of nodes in the approximate k-clique densest subgraph;
  - the number of edges in the approximate k-clique densest subgraph;
  - the edge density of the approximate k-clique densest subgraph;
  - the k-clique density of the approximate k-clique densest subgraph;
  - the computed upper bound on the maximum k-clique density;
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

unsigned CLIQUES_TO_SAMPLE = 10000000;
unsigned sampled_cliques_reserved_size; // Maximum number of cliques for memory allocation; will increase if needed
unsigned *cknodes; // Nodes of a clique being formed
unsigned *ck; // List of all sampled cliques
unsigned *p_ckend; // Pointer to the end of ck[]
unsigned long long cnt_clique; // Number of cliques
unsigned long long cnt_sampled_clique; // Number of sampled cliques
double sampling_prob; // Sampling probability
unsigned long long cnt_clique_in_densest_subgraph; // Number of cliques in the densest subgraph (without sampling)

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
    if (task == SAMPLING)
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
      sampled_cliques_reserved_size = 1.1 * CLIQUES_TO_SAMPLE;
      sampling_prob = (CLIQUES_TO_SAMPLE < cnt_clique) ? (double)CLIQUES_TO_SAMPLE / cnt_clique : 1;
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
      free(cknodes);
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
      // printf("Number of %u-cliques in the densest subgraph: %llu\n", k, cnt_clique_in_densest_subgraph);
      // printf("Density: %.12f\n", (double)cnt_clique_in_densest_subgraph / g->n);
      break;
    }
  }
}

unsigned *perm;
unsigned *rho;
unsigned *ordered_dec_rho;
unsigned *rank_of_rho;
unsigned *rho_pushed_to_max_rank;
bool *is_in_densest_subgraph; // Whether each node is in the densest subgraph

typedef struct {
  unsigned n; // Number of nodes
  unsigned m; // Number of edges
  double density;
  double ub; // An upper bound of maximum density
} DensestSubsetInfo;

static int NodeRhoValueCmp(const void *a, const void *b) {
  return rho[*(const unsigned *)b] - rho[*(const unsigned *)a];
}

void InMemoryFrankWolfe(Graph *g, const unsigned char k) {
  for (int i = cnt_sampled_clique - 1; i >= 0; --i) {
    // Shuffle
    unsigned id = LargeRand() % (i + 1);
    unsigned temp = perm[i];
    perm[i] = perm[id];
    perm[id] = temp;

    // Sequential update
    id = perm[i];
    unsigned node_getting_weight = ck[id * k];
    for (unsigned j = 1; j < k; ++j) {
      if (rho[ck[id * k + j]] < rho[node_getting_weight])
        node_getting_weight = ck[id * k + j];
    }
    ++rho[node_getting_weight];
  }
}

DensestSubsetInfo ExtractDensest(Graph *g, const unsigned char k, unsigned T) {
  // Sort the nodes in decreasing order of rho value
  DensestSubsetInfo info;
  for (unsigned i = 0; i < g->n; ++i) {
    ordered_dec_rho[i] = i;
    rho_pushed_to_max_rank[i] = 0;
  }
  qsort(ordered_dec_rho, g->n, sizeof(unsigned), NodeRhoValueCmp); // Reorder the nodes by decreasing rho values
  for (unsigned i = 0; i < g->n; ++i)
    rank_of_rho[ordered_dec_rho[i]] = i;

  // Iterate over all sampled cliques
  for (unsigned i = 0; i < cnt_sampled_clique; ++i) {
    unsigned node_getting_weight = ck[i * k];
    for (unsigned j = 1; j < k; ++j) {
      if (rank_of_rho[ck[i * k + j]] > rank_of_rho[node_getting_weight])
        node_getting_weight = ck[i * k + j];
    }
    ++rho_pushed_to_max_rank[rank_of_rho[node_getting_weight]];
  }

  // Find the densest subset
  info.density = -1;
  for (unsigned i = 0, cnt_clique_in_subgraph = 0; i < g->n; ++i) {
    cnt_clique_in_subgraph += rho_pushed_to_max_rank[i];
    if (info.density < (double)cnt_clique_in_subgraph / (i + 1)) {
      info.n = i + 1;
      info.m = cnt_clique_in_subgraph;
      info.density = (double)cnt_clique_in_subgraph / (i + 1);
    }
  }
  for (unsigned i = 0; i < info.n; ++i)
    is_in_densest_subgraph[ordered_dec_rho[i]] = true;
  for (unsigned i = info.n; i < g->n; ++i)
    is_in_densest_subgraph[ordered_dec_rho[i]] = false;

  // Compute an upper bound of maximum density
  unsigned sum = 0;
  info.ub = 0;
  double ip1ck = 0; // (i + 1) choose k
  for (unsigned i = 0; i < g->n; ++i) {
    sum += rho[ordered_dec_rho[i]];
    if (i + 1 == k)
      ip1ck = 1;
    else if (i + 1 > k)
      ip1ck = (ip1ck * (i + 1)) / (i + 1 - k);
    if (ip1ck < (double)sum / T)
      info.ub = ip1ck / (i + 1);
    else {
      if (info.ub < (double)sum / T / (i + 1))
        info.ub = (double)sum / T / (i + 1);
      break;
    }
  }
  return info;
}

EdgeList *MakeDensestSubgraphEdgeList(Graph *g, const unsigned char k, const unsigned densest_subset_size) {
  EdgeList *el = (EdgeList *)malloc(sizeof(EdgeList));
  el->n = densest_subset_size;
  el->e = 0;
  for (unsigned i = 0; i < g->e; ++i)
    el->e += (is_in_densest_subgraph[g->edges[i].s] && is_in_densest_subgraph[g->edges[i].t]);
  el->edges = (Edge *)malloc(el->e * sizeof(Edge));
  for (unsigned i = 0, j = 0; i < g->e; ++i) {
    if (is_in_densest_subgraph[g->edges[i].s] && is_in_densest_subgraph[g->edges[i].t]) {
      el->edges[j].s = rank_of_rho[g->edges[i].s];
      el->edges[j].t = rank_of_rho[g->edges[i].t];
      ++j;
    }
  }
  return el;
}

void SampleCliques(Graph *g, const unsigned char k) {
	// Count the number of clqiues
  KCLIST_CliqueEnum(g, k, COUNT);

  // Sampling
  KCLIST_CliqueEnum(g, k, SAMPLING);
}

void Solve(Graph *g, const unsigned char k, unsigned num_iter, clock_t t0) {
  perm = (unsigned *)malloc(cnt_sampled_clique * sizeof(unsigned));
  rho = (unsigned *)calloc(g->n, sizeof(unsigned)); // Initialized to 0 automatically
  is_in_densest_subgraph = (bool *)malloc(g->n * sizeof(bool));
  ordered_dec_rho = (unsigned *)malloc(g->n * sizeof(unsigned));
  rank_of_rho = (unsigned *)malloc(g->n * sizeof(unsigned));
  rho_pushed_to_max_rank = (unsigned *)calloc(g->n, sizeof(unsigned)); // Initialized to 0 automatically
  for (unsigned i = 0; i < cnt_sampled_clique; ++i)
    perm[i] = i;
  for (unsigned T = 1, t = 1; T <= num_iter; T <<= 1) {
    // Step 1: run the Frank-Wolfe based algorithm for num_iter rounds
    for (; t <= T; ++t) {
      if (t % 100 == 0)
        printf("Run round %u...\n", t);
      InMemoryFrankWolfe(g, k);
    }

    // Step 2: give a tentative decomposition
    DensestSubsetInfo info = ExtractDensest(g, k, T);

    // Step 3: count the number of cliques in the densest subset by constructing another Graph
    EdgeList *el = MakeDensestSubgraphEdgeList(g, k, info.n);
    SortByCore(el);
    Relabel(el);
    Graph *p_densest_subgraph = MakeGraph(el);
    KCLIST_CliqueEnum(p_densest_subgraph, k, COUNT_IN_SUBGRAPH);

    //DensestSubsetInfo info = CDF_FindDensestSubset(g, k, T);
    clock_t t1 = clock();
    double edge_density = p_densest_subgraph->e * 2.0 / p_densest_subgraph->n / (p_densest_subgraph->n - 1);
    double density = (double)cnt_clique_in_densest_subgraph / p_densest_subgraph->n;
    double upper_bound = info.ub / sampling_prob / (1 - sqrt(6 * log(g->n) / info.ub));
    printf("Approximate densest subgraph: %u nodes, %u edges, edge density = %f, k-clique density = %f, upper bound = %f. %ld microseconds.\n", p_densest_subgraph->n, p_densest_subgraph->e, edge_density, density, upper_bound, (t1 - t0) * 1000 / CLOCKS_PER_SEC);
    free(p_densest_subgraph);
    //fprintf(ofp, "%u\t%u\t%u\t%.12f\t%.12f\t%ld\n", T, info.n, info.m, info.density, info.ub, t1 - t0);
    fflush(stdout);
  }
  free(perm);
  free(rho);
  free(is_in_densest_subgraph);
  free(ordered_dec_rho);
  free(rank_of_rho);
  free(rho_pushed_to_max_rank);
}

int main(int argc, char **argv) {
  srand(time(NULL));
  EdgeList *el;
  Graph *g;
  unsigned num_threads = atoi(argv[1]);
  unsigned num_iter = atoi(argv[2]);
  unsigned char k = atoi(argv[3]);
  char *file_name = argv[4];
  omp_set_num_threads(num_threads);

  clock_t t0, t1, t2;
  t0 = t1 = clock();

  printf("Reading edgelist from file %s\n", file_name);
  el = ReadEdgeList(file_name);
  printf("Number of nodes = %u\n", el->n);
  printf("Number of edges = %u\n", el->e);
  t2 = clock();
  printf("- Time = %ldh%ldm%lds%ldms\n",(t2 - t1) / CLOCKS_PER_SEC / 3600, ((t2 - t1) / CLOCKS_PER_SEC % 3600) / 60, ((t2 - t1) / CLOCKS_PER_SEC % 60), (t2 - t1) % CLOCKS_PER_SEC * 1000 / CLOCKS_PER_SEC);
  t1 = t2;

  printf("Building the graph structure\n");
  SortByCore(el); // Do core decomposition and render degeneracy ordering to the nodes
  Relabel(el);
  g = MakeGraph(el);
  printf("Number of nodes (degree > 0) = %u\n", g->n);
  t2 = clock();
  printf("- Time = %ldh%ldm%lds%ldms\n", (t2 - t1) / CLOCKS_PER_SEC / 3600, ((t2 - t1) / CLOCKS_PER_SEC % 3600) / 60, ((t2 - t1) / CLOCKS_PER_SEC % 60), (t2 - t1) % CLOCKS_PER_SEC * 1000 / CLOCKS_PER_SEC);
  t1 = t2;

  SampleCliques(g, k);
  t2 = clock();
  printf("- Time = %ldh%ldm%lds%ldms\n", (t2 - t1) / CLOCKS_PER_SEC / 3600, ((t2 - t1) / CLOCKS_PER_SEC % 3600) / 60, ((t2 - t1) / CLOCKS_PER_SEC % 60), (t2 - t1) % CLOCKS_PER_SEC * 1000 / CLOCKS_PER_SEC);
  t1 = t2;

  Solve(g, k, num_iter, t0);
  t2 = clock();
  printf("- Time = %ldh%ldm%lds%ldms\n", (t2 - t1) / CLOCKS_PER_SEC / 3600, ((t2 - t1) / CLOCKS_PER_SEC % 3600) / 60, ((t2 - t1) / CLOCKS_PER_SEC % 60), (t2 - t1) % CLOCKS_PER_SEC * 1000 / CLOCKS_PER_SEC);
  t1 = t2;

  FreeGraph(g);
  printf("- Overall time = %ldh%ldm%lds%ldms\n", (t2 - t0) / CLOCKS_PER_SEC / 3600, ((t2 - t0) / CLOCKS_PER_SEC % 3600) / 60, ((t2 - t0) / CLOCKS_PER_SEC % 60), (t2 - t0) % CLOCKS_PER_SEC * 1000 / CLOCKS_PER_SEC);
  //fprintf(ofp, "%ld\n", t2 - t0);
  //fclose(ofp);
  return 0;
}
